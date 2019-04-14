import sys
import os
import fnmatch
sys.path.append('../tools/')
sys.path.append('../Common/')
from obspy.core import read
from obspy.core.utcdatetime import UTCDateTime
from pyrocko import obspy_compat, io
import Logfile
from ConfigFile import ConfigObj, FilterCfg
import ttt
import obspy
from collections import OrderedDict


def getStation(streamID, MetaDict):

    for i in MetaDict:
        if fnmatch.fnmatch(streamID, i.getName()):
            stobj = i

    return stobj


def getGain(streamID, MetaDict):

    gain = 1

    for i in MetaDict:
        if fnmatch.fnmatch(streamID, i.getName()):
            gain = i.gain

    return gain


def readWaveforms(stationList, tw, EventPath, Origin):

    t2 = UTCDateTime(Origin.time)
    sdspath = os.path.join(EventPath, 'data', str(t2.year))

    Wdict = OrderedDict()

    for i in stationList:
        streamData = i.getName() + '.D.' + str(t2.year)\
                     + '.' + str("%03d" % t2.julday)
        entry = os.path.join(sdspath, i.net, i.sta, i.comp+'.D', streamData)
        tdiff = tw['end']-tw['start']

        try:
            st = read(entry, format="MSEED", starttime=tw['start'],
                      endtime=tw['end'], nearest_sample=True)

        except Exception:
            Logfile.error('readWaveforms: File not found', entry)
            pass

        if len(st.get_gaps()) > 0:
            st.merge(method=0, fill_value='interpolate',
                     interpolation_samples=0)

        if len(st) > 0:
            trdiff = st[0].stats.endtime-st[0].stats.starttime

            totaldiff = abs(trdiff - tdiff)

            if totaldiff < 1:
                Wdict[i.getName()] = st
                Logfile.add(i.getName() + ' added to StreamList ')
        else:
            print(' OUT ', streamData)

    Logfile.red('%d Streams added with available Data' % len(Wdict))
    return Wdict


def readWaveformsPyrocko(stationlist, w, EventPath, Origin, desired):
    Wdict = OrderedDict()
    if desired is 'Z':
        traces = io.load(EventPath+'/data/traces.mseed')
    else:
        traces = io.load(EventPath+'/data/traces_rotated.mseed')
    obspy_compat.plant()

    traces_dict = []
    for il in stationlist:
            for tr in traces:
                tr_name = str(tr.network+'.'+tr.station+'.'+tr.location
                              + '.' + tr.channel[:3])
                if tr_name == str(il) and tr.channel[-1] == desired:
                    st = obspy.Stream()
                    es = obspy_compat.to_obspy_trace(tr)
                    st.extend([es])
                    traces_dict.append(tr)
                    Wdict[il.getName()] = st
    return Wdict


def readWaveformsPyrockodummy(stationlist, w, EventPath, Origin):

    Wdict = OrderedDict()
    for il in stationlist:
        Wdict[il.getName()] = 1.
    return Wdict


def readWaveformsPyrocko_restituted(stationlist, w, EventPath, Origin, desired):
    Wdict = OrderedDict()
    if desired is 'Z':
        traces = io.load(EventPath+'/data/traces_rotated.mseed')
    else:
        traces = io.load(EventPath+'/data/traces_rotated.mseed')

    obspy_compat.plant()

    traces_dict = []
    for tr in traces:
        for il in stationlist:
                tr_name = str(tr.network+'.'+tr.station+'.'+tr.location
                              + '.' + tr.channel[:3])
                if tr_name == str(il) and tr.channel[-1] == desired:
                    st = obspy.Stream()
                    es = obspy_compat.to_obspy_trace(tr)
                    st.extend([es])
                    traces_dict.append(tr)
                    Wdict[il.getName()] = st
                print(tr)
    return Wdict


def readWaveforms_colesseo(stationlist, w, EventPath, Origin, C):
    Wdict = OrderedDict()
    Config = C.parseConfig('config')
    cfg = ConfigObj(dict=Config)
    traces_dict = []
    traces = io.load(cfg.colosseo_scenario_yml()[:-12]+'scenario.mseed')

    for tr in traces:
        for il in stationlist:
                tr_name = str(tr.network+'.'+tr.station+'.'+tr.location
                              + '.' + tr.channel[:3])
                if tr_name == str(il):
                    st = obspy.Stream()
                    es = obspy_compat.to_obspy_trace(tr)
                    st.extend([es])
                    traces_dict.append(tr)
                    Wdict[il.getName()] = st
    return Wdict


def writeWaveform(Folder, station, Stream, flag, network):

    if flag == 'U':
        s1 = 'unfiltered'
    elif flag == 'F':
        s1 = 'filtered'
    elif flag == 'R':
        s1 = 'resampled'

    else:
        Logfile.abort('writeWaveform: Illegal flag <' + flag + '>')

    fname = ('%s-%s-%s.mseed') % (network, station, flag)
    name = os.path.join(Folder['mseed'], fname)
    Stream.write(name, format='MSEED')

    Logfile.add('%s stream for station %s written ' % (s1, station))


def resampleWaveform_2(Waveform, end_frequence):
    return resampleWaveform(Waveform, end_frequence)


def resampleWaveform(Waveform, end_frequence):
    Waveform.resample(end_frequence)
    return Waveform


def resampledummy(Waveform, end_frequence):
    return Waveform


def processWaveforms(WaveformDict, Config, Folder, network, MetaDict, Event,
                     switch, Xcorr):

    Logfile.red('Start Processing')

    cfg = FilterCfg(Config)
    new_frequence = cfg.newFrequency()

    for index, i in enumerate(WaveformDict):
            Logfile.add('%s:%s ---------------------' % (index, i))

            if Config['export_unfiltered'].capitalize() is 'True':
                writeWaveform(Folder, i, WaveformDict[i], 'U', network)

            station = getStation(i, MetaDict)

            if cfg.Int('fm') == 1:
                azi = ttt.bearing(Event, station)
                bazi = ttt.backazi(station, Event)

                msg = 'Takeoff ' + str(station.takeoff) + ' Azi ' + str(azi) +\
                      'Bazi ' + str(bazi)

                Logfile.add(msg)

            gain = float(station.gain)

            if gain == 0.0 or gain == -1.0:
                gain = 1

            WaveformDict[i][0].data * (1.0 / gain)

            if switch is 0:
                Logfile.add('bandpass filtered stream for station %s ' % (i))

                WaveformDict[i].filter('bandpass',
                                       freqmin=cfg.flo(),
                                       freqmax=cfg.fhi(),
                                       corners=cfg.ns(),
                                       zerophase=bool(Config['zph']))

            elif switch is 1:
                Logfile.add('bandpass filtered stream for station %s ' % (i))

                WaveformDict[i].filter('bandpass',
                                       freqmin=cfg.flo2(),
                                       freqmax=cfg.fhi2(),
                                       corners=cfg.ns2(),
                                       zerophase=bool(Config['zph2']))

            else:
                Logfile.add('no filter set for station %s ' % (i))

            if Config['export_filtered'].capitalize() is 'True':
                writeWaveform(Folder, i, WaveformDict[i], 'F', network)

            j = resampleWaveform(WaveformDict[i][0], new_frequence)
            WaveformDict[i] = j

            if Config['export_resampled'].capitalize() is 'True':
                writeWaveform(Folder, i, WaveformDict[i], 'R', network)

    return WaveformDict


def processdummyWaveforms(WaveformDict, Config, Folder, network, MetaDict,
                          Event, switch, Xcorr):
    for index, i in enumerate(WaveformDict):

        WaveformDict[i] = i
    return WaveformDict
