import os
import sys
sys.path.append ('../Common/')

import logging
import sys
import fnmatch

from obspy.core.utcdatetime import UTCDateTime
from obspy.core import read
from obspy.core.stream import Stream,Trace
import obspy.signal.cross_correlation
from obspy.signal.trigger import trigger_onset as triggerOnset
from obspy.signal.trigger import recursive_sta_lta as recSTALTA
from obspy.signal.trigger import classic_sta_lta as classicSTALTA
from obspy.signal.trigger import plot_trigger as plotTrigger
from pyrocko import obspy_compat
from pyrocko import orthodrome, model
from pyrocko import cake
import numpy as np
km = 1000.

import Basic
import Logfile
import Debug
from ObspyFkt   import loc2degrees, obs_TravelTimes
from ConfigFile import ConfigObj, FilterCfg

from config import Trigger                                 # Import from Tools
from waveform import resampleWaveform_2, filterWaveform_2    # Import from Process

logger = logging.getLogger(sys.argv[0])


class Corr (object):

    def __init__(self, shift, value, sampleindex):

        self.shift = shift
        self.value = value
        self.sampleindex = sampleindex


def cmpFilterMetavsXCORR(XcorrMeta, StationMetaList):

    FilterList = []

    for i in StationMetaList:
        for j in XcorrMeta.iterkeys():
            if i.getName() == j:
                FilterList.append(i)

    n1 = len(FilterList)
    n2 = len(StationMetaList)

    Logfile.red('Xcorr Procedure finished %d of %d stations \
    left for processing' % (n1, n2))
    return FilterList


def getArrayShiftValue(refshiftfile, arrayname):

    fobj = open(refshiftfile, 'r')

    for line in fobj:
        line = line.split()

        if fnmatch.fnmatch(line[0], arrayname):
            refshift = line[1]

    fobj.close()
    return refshift


class Xcorr(object):

    def __init__(self, Origin, StationMeta, EventPath, Config, ArrayFolder):

        self.Origin = Origin
        self.StationMeta = StationMeta
        self.EventPath = EventPath
        self.Config = Config
        self.AF = ArrayFolder
        self.mintforerun = int(10)

    # ---------------------------------------------------------------------------------------------

    def calculateTimeWindows(self, mint):

        tw = {}
        st = str(self.Origin.time)[:-1]

        tw['start'] = UTCDateTime(UTCDateTime(st) + (mint - 100))
        tw['end'] = tw['start'] + 200

        tw['xcorrstart'] = UTCDateTime(UTCDateTime(st) + (mint - self.mintforerun))
        tw['xcorrend'] = tw['xcorrstart'] + 20

        Logfile.add(' ORIGIN TIME %s' % UTCDateTime(self.Origin.time))
        Logfile.add(' OVERALL TIME WINDOW : %s - %s' % (tw['start'],      tw['end']))
        Logfile.add(' XCROSS TIME WINDOW  : %s - %s' % (tw['xcorrstart'], tw['xcorrend']))

        return tw

    # ---------------------------------------------------------------------------------------------

    def signoise(self, Waveform, ttime, path):

        st = str(self.Origin.time)[:-1]
        ponset = UTCDateTime(st) + ttime

        winnoise_start = Waveform.stats.starttime+20
        winnoise_end = ponset - 10
        winsig_start = ponset - 2
        winsig_end = ponset + 10

        try:
            winnoise = read(path, format="MSEED", starttime=winnoise_start,
                            endtime=winnoise_end, nearest_sample=True)
            winsig = read(path, format="MSEED", starttime=winsig_start,
                          endtime=winsig_end, nearest_sample=True)
        except:
            Logfile.exception('signoise')

        psignal = abs(winsig.max()[0])
        pnoise = abs(winnoise.max()[0])

        signoise = float(psignal) / float(pnoise)

        return signoise

    # ---------------------------------------------------------------------------------------------

    def resampleWaveform(self, Waveform, end_frequence):

        Logfile.add('enter resampling in crosscorrelation')
        print 'sampling_rate = ', Waveform.stats.sampling_rate
        return resampleWaveform_2(Waveform, end_frequence)

    # ---------------------------------------------------------------------------------------------

    def filterWaveform(self, Waveform):

        Logfile.red('Filter Waveform: ')
        cfg = FilterCfg(self.Config)

        new_frequence = (cfg.newFrequency())

        st = Stream()
        for i in Waveform:
                Logfile.red('Downsampling to %s: from %d' % (new_frequence,
                            i.stats.sampling_rate))

                j = i.resample(new_frequence)

                switch = cfg.filterswitch()

                if switch == 1:
                        Logfile.add('bandpass filtered \
                        stream for station %s ' % (i))

                        j.filter('bandpass',
                                 freqmin=cfg.flo(),
                                 freqmax=cfg.fhi(),
                                 corners=cfg.ns(),
                                 zerophase=bool(self.Config['zph']))

                elif switch == 2:
                        Logfile.add('bandpass filtered \
                        stream for station %s ' % (i))

                        j.filter('bandpass',
                                 freqmin=cfg.flo2(),
                                 freqmax=cfg.fhi2(),
                                 corners=cfg.ns2(),
                                 zerophase=bool(self.Config['zph']))
                st.append(j)

        return st

    def filterpyrockoWaveform(self, Waveform):
        obspy_compat.plant()
        Logfile.red('Filter Waveform: ')
        cfg = FilterCfg(self.Config)

        new_frequence = (cfg.newFrequency())
        st = Stream()
        for i in Waveform:
            tr = obspy_compat.to_pyrocko_trace(i)
            tr.downsample_to(new_frequence)

            Logfile.red('Downsampling to %d: from %d' % (new_frequence,
                        i.stats.sampling_rate))

            if switch == 0:
                Logfile.add('bandpass filtered stream for station %s ' % (i))

                tr.bandpass(4, cfg.flo, cfg.fhi)
            elif switch == 1:
                Logfile.add('bandpass filtered stream for station %s ' % (i))
            j = obspy_compat.to_obspy_trace(tr)

            st.append(j)

        return st

    # ---------------------------------------------------------------------------------------------

    def readWaveformsCross(self, station, tw, ttime):

        t2 = UTCDateTime(self.Origin.time)
        sdspath = os.path.join(self.EventPath, 'data', str(t2.year))

        stream = ''
        snr = ''

        if station.loc == '--':
            station.loc = ''

        streamData = station.net + '.' + station.sta + '.' + station.loc + '.' + station.comp + '.D.' + str(t2.year) + '.' + str("%03d" % t2.julday)
        entry = os.path.join(sdspath, station.net, station.sta, station.comp + '.D', streamData)
        st= read(entry, format="MSEED", starttime=tw['start'], endtime=tw['end'], nearest_sample=True)
        if len(st.get_gaps()) > 0:
            st.merge (method=0, fill_value='interpolate', interpolation_samples=0)
        snr  = self.signoise     (st[0], ttime, entry)
        #amp = self.maxAmplitude (st[0], ttime, Origin, entry)
        try:
            stream = self.filterWaveform(st)
        except:
            stream = self.filterWaveform(st)
    #        stream = self.filterpyrockoWaveform(st)

       # stream = st
        xname  = os.path.join(self.AF,(streamData+'_all.mseed'))
        stream.write (xname,format='MSEED')
        stream.trim (tw['xcorrstart'], tw['xcorrend'])
      #  stream[0].stats.starttime = UTCDateTime(3600)
        #stream.write(streamData,format='MSEED')
      #  stream[0].stats.starttime = UTCDateTime(1971, 1, 1, 1, 0)
        return stream, snr

    # ---------------------------------------------------------------------------------------------

    def traveltimes(self):

        Logfile.red ('Enter AUTOMATIC CROSSCORRELATION ')
        Logfile.red ('\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++\n ')
        T     = []
        Wdict = {}
        SNR   = {}

        for i in self.StationMeta:

            Logfile.red ('read in %s '%(i))
            de = loc2degrees     (self.Origin, i)
    	    Phase = cake.PhaseDef('P')
            model = cake.load_model()
            arrivals= model.arrivals([de,de], phases=Phase, zstart=self.Origin.depth*km)
	    try:
                	ptime = arrivals[0].t
	    except:
			try:
		        	arrivals= model.arrivals([de,de], phases=Phase, zstart=o_depth*km-2.1)
				ptime = arrivals[0].t
			except:
				ptime = ptime
	    T.append(ptime)
            if ptime == 0:
                Logfile.red ('Available phases for station %s in range %f deegree' % (i,de))
                Logfile.red ('you tried phase %s' % (self.Config[phasename]))
                raise Exception ("ILLEGAL: phase definition")

            tw = self.calculateTimeWindows(ptime)
            try:
                w, snr = self.readWaveformsCross (i, tw, ptime)
                Wdict [i.getName()] = w
                SNR   [i.getName()] = snr
            except:
                pass

            Logfile.red ('\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++ ')

        Logfile.red ('Exit AUTOMATIC FILTER ')
        print SNR
        return Wdict, SNR

    # ---------------------------------------------------------------------------------------------

    def readWaveformsPicker (self,station, tw, Origin, ttime):

        t2      = UTCDateTime(self.Origin.time)
        sdspath = os.path.join(self.EventPath,'data', str(t2.year))

        if station.loc == '--':
            station.loc = ''

        staName    = station.net + '.'   + station.sta  + '.' + station.loc + '.' + station.comp
        streamData = staName     + '.D.' + str(t2.year) + '.' + str("%03d" % t2.julday)
        entry      = os.path.join (sdspath, station.net, station.sta, station.comp + '.D', streamData)
        st         = read (entry, format="MSEED", starttime=tw['start'], endtime=tw['end'], nearest_sample=True)

        if len(st.get_gaps()) > 0:
            st.merge (method=0, fill_value='interpolate', interpolation_samples=0)

        stream = self.filterWaveform (st)
        return stream

    # ---------------------------------------------------------------------------------------------

    def searchMeta (self,sname,Metalist):

        for i in Metalist:
            if sname == i.getName():
                return i
    # ---------------------------------------------------------------------------------------------

    def refTrigger (self,RefWaveform):

        name = ('%s.%s.%s.%s') % (RefWaveform[0].stats.network, RefWaveform[0].stats.station,
                                  RefWaveform[0].stats.location,RefWaveform[0].stats.channel)

        i     = self.searchMeta (name,self.StationMeta)
        de    = loc2degrees     (self.Origin, i)
        tt    = obs_TravelTimes (de, self.Origin.depth)
        ptime = 0

        Phase = cake.PhaseDef('P')
        model = cake.load_model()

        arrivals= model.arrivals([de,de], phases=Phase, zstart=self.Origin.depth*km)
        try:
            ptime = arrivals[0].t
        except:
            arrivals= model.arrivals([de,de], phases=Phase, zstart=o_depth*km-0.1)
            ptime = arrivals[0].t
        phasename = ('%sphase') % (os.path.basename(self.AF))

        if ptime == 0:
                print '\033[31mAvailable phases for reference station %s in range %f deegree\033[0m'%(i,de)
                print '\033[31m'+'|'.join([str(item['phase_name']) for item in tt])+'\033[0m'
                print '\033[31myou tried phase %s\033[0m'%(self.Config[phasename])
                raise Exception("\033[31mILLEGAL: phase definition\033[0m")

        tw  = self.calculateTimeWindows (ptime)
        stP = self.readWaveformsPicker  (i, tw, self.Origin, ptime)

        refuntouchname = os.path.basename(self.AF)+'-refstation-raw.mseed'
        stP.write (os.path.join(self.EventPath,refuntouchname),format='MSEED',byteorder='>')
        print float(self.Config['refstationfreqmin'])
        print float(self.Config['refstationfreqmax'])
        stP.filter("bandpass", freqmin  = float (self.Config['refstationfreqmin']),
                               freqmax  = float (self.Config['refstationfreqmax']))

        stP.trim(tw['xcorrstart'], tw['xcorrend'])
        trP = stP[0]

        trP.stats.starttime = UTCDateTime(3600)
        refname = os.path.basename (self.AF)+'-refstation-filtered.mseed'
        trP.write (os.path.join (self.EventPath,refname),format='MSEED',byteorder='>')

        sta = float (self.Config['refsta'])
        lta = float (self.Config['reflta'])
        cft = recSTALTA (trP.data, int(sta * trP.stats.sampling_rate), int (lta * trP.stats.sampling_rate))

        t = triggerOnset (cft,lta,sta)

        try:
            onset = t[0][0] / trP.stats.sampling_rate
            print 'ONSET ',onset

        except:
            onset = self.mintforerun

        trigger = trP.stats.starttime+onset

        print 'TRIGGER ',trigger
        print 'THEORETICAL: ',UTCDateTime(3600)+self.mintforerun
        tdiff = (trP.stats.starttime+onset)-(UTCDateTime(3600)+self.mintforerun)
        print 'TDIFF: ',tdiff

        refp = UTCDateTime(self.Origin.time)+ptime
        reftriggeronset = refp+onset-self.mintforerun

        if int(self.Config['autoxcorrcorrectur']) == 1:
            try:

                refmarkername     = os.path.join (self.EventPath,('%s-marker') % (os.path.basename(self.AF)))
                fobjrefmarkername = open(refmarkername,'w')
                fobjrefmarkername.write ('# Snuffler Markers File Version 0.2\n')
                fobjrefmarkername.write (('phase: %s 0 %s    None           None         None         XWStart        None False\n')%(tw['xcorrstart'].strftime('%Y-%m-%d %H:%M:%S.%f'),name))
                fobjrefmarkername.write (('phase: %s 0 %s    None           None         None         XWEnd        None False\n')%(tw['xcorrend'].strftime('%Y-%m-%d %H:%M:%S.%f'),name))
                fobjrefmarkername.write (('phase: %s 1 %s    None           None         None         TheoP        None False\n')%(refp.strftime('%Y-%m-%d %H:%M:%S.%f'),name))
                fobjrefmarkername.write (('phase: %s 3 %s    None           None         None         XTrig        None False')%(reftriggeronset.strftime('%Y-%m-%d %H:%M:%S.%f'),name))
                fobjrefmarkername.close ()

                cmd = 'snuffler %s --markers=%s&'%(os.path.join(self.EventPath,refuntouchname),refmarkername)
                os.system(cmd)

                thrOn  = float (self.Config['reflta']) # 4
                thrOff = float (self.Config['refsta']) # 0.7
                plotTrigger (trP, cft, thrOn, thrOff)


                selection = float (raw_input('Enter self picked phase in seconds: '))
                tdiff     = selection-self.mintforerun
                print selection-self.mintforerun

                refname = os.path.basename (self.AF)+'-shift.mseed'
                trP.stats.starttime = trP.stats.starttime - selection
                trP.write (os.path.join (self.EventPath,refname),format='MSEED')

            except:
                selection = 0.
                refname   = os.path.basename (self.AF)+'-shift.mseed'
                trP.stats.starttime = trP.stats.starttime - selection -self.mintforerun
                trP.write (os.path.join (self.EventPath,refname),format='MSEED')
        '''
        tdiff = 0
        trigger = trP.stats.starttime
        '''
        To = Trigger(name,trigger,os.path.basename(self.AF),tdiff)


        return tdiff,To

    # ---------------------------------------------------------------------------------------------

    def shiftSeismograms (self,StreamDict, XcorrDict,pickerShift):

        L = []
        S = []

        dsfactor = float (self.Config ['xcorrtreshold'])

        for stream in StreamDict.iterkeys():
            for shift in XcorrDict.iterkeys():
                if stream == shift:

                    StreamDict[stream][0].stats.starttime = StreamDict[stream][0].stats.starttime + XcorrDict[shift].shift+pickerShift

                    if XcorrDict[shift].value > dsfactor:
                        fname = stream + '.new'
                        StreamDict[stream].write(os.path.join(self.AF,fname), format='MSEED')

                        if XcorrDict[shift].value < 0:
                            t = -1
                        else:
                            t = 1

                        info = [stream, XcorrDict[shift].shift,t]

                        L.append (info)
                        S.append (stream)

                    if abs(XcorrDict[shift].value) < dsfactor:
                        print 'OUT: ', stream, XcorrDict[shift].value, XcorrDict[shift].sampleindex, XcorrDict[shift].shift
        return L, S

    # ---------------------------------------------------------------------------------------------

    def writeShift (self,ShiftList):
        import csv
        filename = os.path.join (self.AF,'shift.dat')
        #np.savetxt(filename, ShiftList)
        with open(filename, 'wb') as csv_file:
            writer = csv.writer(csv_file)
            for val in ShiftList:
                writer.writerow([val])

        #fobj = open (os.path.join (self.AF,'shift.dat'), 'w')
        #fobj.write ('\n'.join(ShiftList))
        #fobj.close ()

    # ---------------------------------------------------------------------------------------------

    def f6 (self,d1):
        return max (d1, key=d1.get)

    # ---------------------------------------------------------------------------------------------

    def filterCorrDict(self,CorrDict,onset):

        fCD = {}

        dsfactor = float(self.Config['xcorrtreshold'])

        for stream in CorrDict.iterkeys():
            if CorrDict[stream].value >= dsfactor:
                fCD[stream] = CorrDict[stream]
                fCD[stream].value = fCD[stream].value

        return fCD

    # ---------------------------------------------------------------------------------------------

    def doXcorr (self):
        StreamDict, SNRDict = self.traveltimes()
        t = self.f6(SNRDict)
        Logfile.add ('doXcorr: REFERENCE: ' + t)

        for i in SNRDict.iterkeys():
            Logfile.add ('doXcorr: STREAM: ' + i + ' SNR: ' + str(SNRDict[i]))

       #alternativeref = os.path.join (*self.AF.split('/')[-1:])    + 'refstation'   #hs
        alternativeref = os.path.join (*self.AF.split(os.sep)[-1:]) + 'refstation'   #hs

        if self.Config [alternativeref] == '' : t = t
        else:                                   t = self.Config [alternativeref]

        corrDict = {}
        ref      = StreamDict[t][0].data

        Logfile.red ('Reference Station of %s for Xcorr Procedure %s' % (os.path.basename(self.AF),t))
        Logfile.red ('Enter Xcorr Procedure ')

        for stream in StreamDict.iterkeys():
           xcorrshiftvalue = StreamDict[stream][0].stats.npts/10

           a, b  = obspy.signal.cross_correlation.xcorr (ref, StreamDict[stream][0], xcorrshiftvalue)
           shift = a / StreamDict[stream][0].stats.sampling_rate
           corrDict[stream] = Corr(shift, b, a)

           msg =  'Index: ' + str(a) + ' Value: ' +str(b) + ' ----> '
           msg += (str(stream) + str(StreamDict[stream][0].stats.sampling_rate) + ' SHIFT IN TIME: ' + str(shift))
           Logfile.add (msg)

        Logfile.red ('Finish Xcorr Procedure ')
        return corrDict,StreamDict[t],StreamDict

    # ---------------------------------------------------------------------------------------------

    def runXcorr (self):

        CD,ref,WD = self.doXcorr()
        onset     = 0
        #onset    = self.refTrigger (ref,self.EventPath,self.StationMeta)
        tdiff,triggerobject = self.refTrigger(ref)

        fCD = self.filterCorrDict(CD,onset)

        C, Streams = self.shiftSeismograms (WD, fCD,onset)

        self.writeShift(C)
        return fCD,triggerobject
