import sys
import logging
import os.path as op
from optparse import OptionParser

from pyrocko import util, scenario, guts, gf

import os
import sys
import logging
import fnmatch
sys.path.append('../tools/')
sys.path.append('../Common/')
from obspy.core             import read
from obspy.core.stream      import Stream, Trace
from obspy.core.utcdatetime import UTCDateTime
from pyrocko import obspy_compat, io
from pyrocko import orthodrome, model
import numpy
import  Basic
import  Globals
import  Logfile
from    ConfigFile import ConfigObj, FilterCfg, SynthCfg
from    config import Station,Event     #       Import from Tools
import  ttt                             #       Import from Process
import obspy
from collections import OrderedDict

# -------------------------------------------------------------------------------------------------

def getStation(streamID,MetaDict):

    for i in MetaDict:
        if fnmatch.fnmatch(streamID,i.getName()) :  stobj = i

    return stobj

# -------------------------------------------------------------------------------------------------

def getGain(streamID,MetaDict):

    gain=1

    for i in MetaDict:
        if fnmatch.fnmatch(streamID,i.getName()):
            gain = i.gain

    return gain
# -------------------------------------------------------------------------------------------------

def readWaveforms(stationList,tw,EventPath,Origin):

    t2  = UTCDateTime(Origin.time)
    sdspath = os.path.join(EventPath,'data',str(t2.year))

    Wdict = OrderedDict()

    for i in stationList:
        streamData = i.getName() + '.D.' + str(t2.year) + '.' + str("%03d" % t2.julday)
        entry  = os.path.join(sdspath,i.net,i.sta,i.comp+'.D',streamData)
        tdiff  = tw['end']-tw['start']

        try:
            st = read(entry,format="MSEED",starttime=tw['start'],endtime=tw['end'],nearest_sample=True)

        except:
            Logfile.error('readWaveforms : File not found', entry)
            pass

        if len(st.get_gaps()) > 0:
            st.merge(method=0, fill_value='interpolate', interpolation_samples=0)


        if len(st) > 0:
            trdiff = st[0].stats.endtime-st[0].stats.starttime

            totaldiff = abs(trdiff - tdiff)

            if totaldiff < 1:
                Wdict[i.getName()] = st
                Logfile.add(i.getName() + ' added to StreamList ')
        else:
           print(' OUT ',streamData)

    Logfile.red('%d Streams added with available Data' % len(Wdict))
    return Wdict

def readWaveformsPyrocko(stationlist, w,EventPath,Origin):
    Wdict = OrderedDict()
    traces = io.load(EventPath+'/data/traces.mseed')
    obspy_compat.plant()

    traces_dict = []
    for il in stationlist:
            for tr in traces:
              tr_name = str(tr.network+'.'+tr.station+'.'+tr.location+'.'+tr.channel[:3])
              if tr_name == str(il):
                    st = obspy.Stream()
                    es = obspy_compat.to_obspy_trace(tr)
                    st.extend([es])
                    traces_dict.append(tr)
                    Wdict[il.getName()] = st
    return Wdict

def readWaveformsPyrockodummy(stationlist, w,EventPath,Origin):

    Wdict = OrderedDict()
    for il in stationlist:
        Wdict[il.getName()] = 1.
    return Wdict

def readWaveformsPyrocko_restituted(stationlist, w, EventPath, Origin):
    Wdict = OrderedDict()
    traces = io.load(EventPath+'/data/traces_restituted.mseed')
    obspy_compat.plant()

    traces_dict = []
    for tr in traces:
        for il in stationlist:
                tr_name = str(tr.network+'.'+tr.station+'.'+tr.location+'.'+tr.channel[:3])
                if tr_name == str(il):
                    st = obspy.Stream()
                    es = obspy_compat.to_obspy_trace(tr)
                    st.extend([es])
                    traces_dict.append(tr)
                    Wdict[il.getName()] = st
    return Wdict


def readWaveforms_colesseo(stationlist, w, EventPath, Origin, C):
    Wdict = OrderedDict()
    Config = C.parseConfig('config')
    cfg = ConfigObj(dict=Config)
    traces_dict = []
    traces = io.load(cfg.colosseo_scenario_yml()[:-12]+'scenario.mseed')

    for tr in traces:
        for il in stationlist:
                tr_name = str(tr.network+'.'+tr.station+'.'+tr.location+'.'+tr.channel[:3])
                if tr_name == str(il):
                    st = obspy.Stream()
                    es = obspy_compat.to_obspy_trace(tr)
                    st.extend([es])
                    traces_dict.append(tr)
                    Wdict[il.getName()] = st
    return Wdict


# -------------------------------------------------------------------------------------------------

def writeWaveform(Folder, station, Stream, flag, network):

    if   flag == 'U':  s1 = 'unfiltered'
    elif flag == 'F':  s1 = 'filtered'
    elif flag == 'R':  s1 = 'resampled'

    else : Logfile.abort('writeWaveform : Illegal flag <' + flag +'>')

    fname =('%s-%s-%s.mseed') %(network,station,flag)
    name  = os.path.join(Folder['mseed'],fname)
    Stream.write(name, format='MSEED')

    Logfile.add ('%s stream for station %s written '%(s1, station))

# -------------------------------------------------------------------------------------------------

def resampleWaveform_2(Waveform, end_frequence):
    return resampleWaveform(Waveform, end_frequence)

def resampleWaveform(Waveform, end_frequence):

          Waveform.resample(end_frequence)
          return Waveform

def resampledummy(Waveform, end_frequence):

          return Waveform
# -------------------------------------------------------------------------------------------------

def filterWaveform_2(Config, wafeform, station) :
    return filterWaveform(Config, wafeform, station)

def filterWaveform(Config, waveform, station) :

    cfg= FilterCfg(self.Config)
    switch = cfg.filterswitch()              # 'filterswitch'

    if switch == 1:
       which = 'bandpass'
       waveform.filter(which, freqmin   = cfg.flo(),           # ['flo']
                               freqmax   = cfg.fhi(),           # ['fhi']
                               corners   = cfg.ns(),            # ['ns']
                               zerophase = bool(Config['zph']))
    elif switch == 2:
       waveform.filter(which, freqmin   = cfg.flo2(),           # ['flo']
                               freqmax   = cfg.fhi2(),           # ['fhi']
                               corners   = cfg.ns2(),            # ['ns']
                               zerophase = bool(Config['zph2']))
    elif switch == 3:
      which = 'highpass'
      waveform.filter(which, freq = cfg.h_fc(),             # ['h_fc']
                              corners  = cfg.h_ns(),             # ['h_ns']
                              zerophase= bool(Config['h_zph']))
    else : return None

    Logfile.add('%s filtered stream for station %s '%(which, i))
    return waveform

# -------------------------------------------------------------------------------------------------

def processWaveforms(WaveformDict,Config,Folder,network,MetaDict,Event,switch,Xcorr):

    Logfile.red('Start Processing')

    cfg  = FilterCfg(Config)
    new_frequence = cfg.newFrequency()                #  ['new_frequence']

    vp,vs,rho = ttt.calcak135parameter(Event)

    for index,i in enumerate(WaveformDict):
            Logfile.add('%s:%s -------------------------------------------------------------' %(index,i))

            if Config ['export_unfiltered'].capitalize() == 'True':
               writeWaveform(Folder,i,WaveformDict[i],'U',network)

            station = getStation(i,MetaDict)

            psign = 1


            if cfg.Int('fm') == 1:
                azi   = ttt.bearing(Event, station)
                bazi  = ttt.backazi(station, Event)
                psign = ttt.dubcup (rho, vp,vs,Event.strike,Event.dip,Event.rake, azi, station.takeoff)

                msg = 'Takeoff ' + str(station.takeoff) + ' Azi ' + str(azi) + ' Bazi ' + str(bazi)
                msg +=(' vp ' + str(vp) + ' vs ' + str(vs) + ' rho ' + str(rho) + ' Psign ' + str(psign))
                Logfile.add(msg)

            try:
            	Logfile.add(station.getName() + ' ' + station.lat + ' ' + station.lon + ' ' + station.gain + ' PSIGN: ' + str(psign))
            except:
                pass

            if psign == -1:
                Logfile.add('correcting polarisation for station %s ' %(i))
                -1.*WaveformDict[i][0].data
            #gain correctur
            gain = float(station.gain)

            if gain == 0.0 or gain == -1.0 : gain = 1    #hs : gain missing in metafile

            WaveformDict[i][0].data *(1.0 / gain)

            if switch == 0:
                Logfile.add('bandpass filtered stream for station %s '%(i))

                WaveformDict[i].filter('bandpass',
                                        freqmin   = cfg.flo(),           # ['flo']
                                        freqmax   = cfg.fhi(),           # ['fhi']
                                        corners   = cfg.ns(),            # ['ns']
                                        zerophase = bool(Config['zph']))

            elif switch == 1:
                Logfile.add('bandpass filtered stream for station %s '%(i))

                WaveformDict[i].filter('bandpass',
                                        freqmin   = cfg.flo2(),           # ['flo']
                                        freqmax   = cfg.fhi2(),           # ['fhi']
                                        corners   = cfg.ns2(),            # ['ns']
                                        zerophase = bool(Config['zph2']))

            else :
               Logfile.add('no filter set for station %s '%(i))

            if Config ['export_filtered'].capitalize() == 'True':
               writeWaveform(Folder,i,WaveformDict[i],'F',network)

            #resampling
            j = resampleWaveform(WaveformDict[i][0],new_frequence)
            WaveformDict[i] = j

            if Config ['export_resampled'].capitalize() == 'True':
               writeWaveform(Folder,i,WaveformDict[i],'R',network)

    return WaveformDict


def processdummyWaveforms(WaveformDict,Config,Folder,network,MetaDict,Event,switch,Xcorr):


    for index,i in enumerate(WaveformDict):

        #resampling
        #j = resampledummy(WaveformDict[i][0],2.)
        WaveformDict[i] = i

    return WaveformDict


def processpyrockoWaveforms(WaveformDict,Config,Folder,network,MetaDict,Event,switch,Xcorr):
    WaveformDict_obs = []
    obspy_compat.plant()
    Logfile.red('Start Processing')
    cfg  = FilterCfg(Config)
    new_frequence = cfg.newFrequency()                #  ['new_frequence']

    traces = []
    for tr in WaveformDict:
            Logfile.add('%s:%s -------------------------------------------------------------' %(index,i))

            if Config ['export_unfiltered'].capitalize() == 'True':
                tr = obspy_compat.to_obspy_trace(tr)
                writeWaveform(Folder,i,tr,'U',network)


            if switch == 0:
                Logfile.add('bandpass filtered stream for station %s '%(i))

                tr.bandpass(4, cfg.flo,cfg.fhi)
            elif switch == 1:
                Logfile.add('bandpass filtered stream for station %s '%(i))

                tr.bandpass(4, cfg.flo2,cfg.fhi2)
            tr.downsample_to(new_frequence)
            tr = obspy_compat.to_obspy_trace(tr)
            traces.append(tr)
            stream = obspy.Stream()
            Wdict= stream.extend([tr])
            WaveformDict_obs.append(Wdict)

    return WaveformDict_obs
