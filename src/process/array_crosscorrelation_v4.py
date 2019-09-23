import sys
import logging
import os
import fnmatch
from obspy.core.utcdatetime import UTCDateTime
from obspy.core import read
from obspy.core.stream import Stream
import obspy.signal.cross_correlation
from obspy.signal.trigger import trigger_onset as triggerOnset
from obspy.signal.trigger import recursive_sta_lta as recSTALTA
from obspy.signal.trigger import plot_trigger as plotTrigger
from pyrocko import obspy_compat,  cake, io
import numpy as num
km = 1000.
from palantiri.common import Logfile
from palantiri.common.ObspyFkt import loc2degrees
from palantiri.common.ConfigFile import ConfigObj, FilterCfg
from palantiri.tools.config import Trigger
from palantiri.process.waveform import resampleWaveform_2
from collections import OrderedDict

logger = logging.getLogger(sys.argv[0])


class Corr(object):

    def __init__(self, shift, value, sampleindex):

        self.shift = shift
        self.value = value
        self.sampleindex = sampleindex


def cmpFilterMetavsXCORR(XcorrMeta, StationMetaList):

    FilterList = []

    for i in StationMetaList:
        if sys.version_info.major >= 3:
            for j in sorted(XcorrMeta.keys()):
                if i.getName() == j:
                    FilterList.append(i)
        else:
            for j in XcorrMeta.keys():
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

    def __init__(self, Origin, StationMeta, EventPath, Config, Syn_in,
                 ArrayFolder):

        self.Origin = Origin
        self.StationMeta = StationMeta
        self.EventPath = EventPath
        self.Config = Config
        self.AF = ArrayFolder
        self.mintforerun = int(10)
        self.Syn_in = Syn_in

    def calculateTimeWindows(self, mint):

        tw = {}
        st = str(self.Origin.time)[:-1]

        tw['start'] = UTCDateTime(UTCDateTime(st) + (mint - 100))
        tw['end'] = tw['start'] + 200

        tw['xcorrstart'] = UTCDateTime(UTCDateTime(st)
                                       + (mint - self.mintforerun))
        tw['xcorrend'] = tw['xcorrstart'] + 40

        Logfile.add(' ORIGIN TIME %s' % UTCDateTime(self.Origin.time))
        Logfile.add(' OVERALL TIME WINDOW: %s - %s' % (tw['start'], tw['end']))
        Logfile.add(' XCROSS TIME WINDOW: %s - %s' % (tw['xcorrstart'],
                                                      tw['xcorrend']))

        return tw

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
        except Exception:
            Logfile.exception('signoise')

        psignal = abs(winsig.max()[0])
        pnoise = abs(winnoise.max()[0])

        signoise = float(psignal) / float(pnoise)

        return signoise

    def resampleWaveform(self, Waveform, end_frequence):

        Logfile.add('enter resampling in crosscorrelation')
        print('sampling_rate = ', Waveform.stats.sampling_rate)
        return resampleWaveform_2(Waveform, end_frequence)

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

    def readWaveformsCross(self, station, tw, ttime):

        t2 = UTCDateTime(self.Origin.time)
        sdspath = os.path.join(self.EventPath, 'data', str(t2.year))

        stream = ''
        snr = ''
        if station.loc == '--':
            station.loc = ''

        streamData = station.net + '.' + station.sta + '.' + station.loc + '.'\
                                 + station.comp + '.D.' + str(t2.year) + '.'\
                                 + str("%03d" % t2.julday)
        entry = os.path.join(sdspath, station.net, station.sta, station.comp +
                             '.D', streamData)
        st = read(entry, format="MSEED", starttime=tw['start'],
                  endtime=tw['end'], nearest_sample=True)
        if len(st.get_gaps()) > 0:
            st.merge(method=0, fill_value='interpolate',
                     interpolation_samples=0)
        snr = self.signoise(st[0], ttime, entry)
        stream = self.filterWaveform(st)

        xname = os.path.join(self.AF, (streamData+'_all.mseed'))
        stream.write(xname, format='MSEED')
        stream.trim(tw['xcorrstart'], tw['xcorrend'])

        return stream, snr

    def readWaveformsCross_pyrocko(self, station, tw, ttime, traces):
        obspy_compat.plant()

        cfg = ConfigObj(dict=self.Config)

        t2 = UTCDateTime(self.Origin.time)

        found = False

        for tr in traces:
            tr_name = str(tr.network+'.'+tr.station+'.'+tr.location +
                          '.'+tr.channel[:3])
            if tr_name == str(station)[:-2] or tr_name == str(station)[:]:
                traces_station = tr
                es = obspy_compat.to_obspy_trace(traces_station)
                streamData = station.net + '.' + station.sta + '.'\
                                         + station.loc + '.'\
                                         + station.comp\
                                         + '.D.'\
                                         + str(t2.year) + '.'\
                                         + str("%03d" % t2.julday)

                st = obspy.Stream()
                st.extend([es])
                stream = ''
                snr = ''
                if station.loc == '--':
                    station.loc = ''

                if len(st.get_gaps()) > 0:
                    st.merge(method=0, fill_value='interpolate',
                             interpolation_samples=0)
                snr_trace = traces_station.chop(tmin=traces_station.tmin,
                                                tmax=traces_station.tmin
                                                + ttime-20.,
                                                inplace=False)
                snr = num.var(snr_trace.ydata)
                stream = self.filterWaveform(st)

                xname = os.path.join(self.AF, (streamData+'_all.mseed'))
                stream.write(xname, format='MSEED')
                stream.trim(tw['xcorrstart'], tw['xcorrend'])
                found = True
                return stream, snr, found
        if found is False:
                    print('Waveform missing!', tr_name, str(station))

    def readWaveformsCross_colesseo(self, station, tw, ttime):
        obspy_compat.plant()
        Config = self.Config
        cfg = ConfigObj(dict=Config)
        t2 = UTCDateTime(self.Origin.time)

        traces = io.load(cfg.colosseo_scenario_yml()[:-12]+'scenario.mseed')

        for tr in traces:
            tr_name = str(tr.network+'.'+tr.station+'.'+tr.location+'.'
                                    + tr.channel[:3])
            if tr_name == str(station):
                    traces_station = tr

                    es = obspy_compat.to_obspy_trace(traces_station)
                    streamData = station.net + '.' + station.sta + '.'\
                                             + station.loc + '.'\
                                             + station.comp + '.D.'\
                                             + str(t2.year) + '.'\
                                             + str("%03d" % t2.julday)

                    st = obspy.Stream()
                    st.extend([es])
                    stream = ''
                    snr = ''

                    if station.loc == '--':
                        station.loc = ''

                    if len(st.get_gaps()) > 0:
                        st.merge(method=0, fill_value='interpolate',
                                 interpolation_samples=0)
                    snr_trace = traces_station.chop(tmin=traces_station.tmin,
                                                    tmax=traces_station.tmin +
                                                    ttime-20.,
                                                    inplace=False)
                    snr = num.var(snr_trace.ydata)
                    stream = self.filterWaveform(st)

                    xname = os.path.join(self.AF, (streamData+'_all.mseed'))
                    stream.write(xname, format='MSEED')
                    stream.trim(tw['xcorrstart'], tw['xcorrend'])
                    return stream, snr

            else:
                pass

    def traveltimes(self, phase, traces):

        Logfile.red('Enter AUTOMATIC CROSSCORRELATION ')
        Logfile.red('\n\n+++++++++++++++++++++++++++++++++++++++++++++++\n ')
        T = []
        Wdict = OrderedDict()
        SNR = OrderedDict()
        Config = self.Config
        cfg = ConfigObj(dict=Config)

        for i in self.StationMeta:
            Logfile.red('read in %s ' % (i))
            de = loc2degrees(self.Origin, i)
            Phase = cake.PhaseDef(phase)
            traveltime_model = cfg.Str('traveltime_model')
            import palantiri
            path = palantiri.__path__
            model = cake.load_model(path[0]+'/data/'+traveltime_model)
            if cfg.colesseo_input() is True:
                arrivals = model.arrivals([de, de], phases=Phase,
                                          zstart=self.Origin.depth, zstop=0.)
            else:
                arrivals = model.arrivals([de, de], phases=Phase,
                                          zstart=self.Origin.depth*km,
                                          zstop=0.)
            try:
                ptime = arrivals[0].t
            except Exception:
                try:
                    arrivals = model.arrivals([de, de], phases=Phase,
                                              zstart=self.Origin.depth*km-2.1)
                    ptime = arrivals[0].t
                except Exception:
                    ptime = 0
            T.append(ptime)
            if ptime == 0:
                Logfile.red('Available phases for station %s in\
                            range %f deegree' % (i, de))
                Logfile.red('you tried phase %s' % (phase))
                raise Exception("ILLEGAL: phase definition")
            else:
                tw = self.calculateTimeWindows(ptime)
                if cfg.pyrocko_download() is True:
                    w, snr, found = self.readWaveformsCross_pyrocko(i, tw, ptime,
                                                             traces)
                elif cfg.colesseo_input() is True:
                    w, snr = self.readWaveformsCross_colesseo(i, tw, ptime)
                else:
                    w, snr = self.readWaveformsCross(i, tw, ptime)
                Wdict[i.getName()] = w
                SNR[i.getName()] = snr

            Logfile.red('\n\n+++++++++++++++++++++++++++++++++++++++++++++++ ')

        Logfile.red('Exit AUTOMATIC FILTER ')
        return Wdict, SNR

    def readWaveformsPicker(self, station, tw, Origin, ttime):

        t2 = UTCDateTime(self.Origin.time)
        sdspath = os.path.join(self.EventPath, 'data', str(t2.year))

        if station.loc == '--':
            station.loc = ''

        staName = station.net + '.' + station.sta + '.' + station.loc\
                              + '.' + station.comp
        streamData = staName + '.D.' + str(t2.year) + '.'\
                             + str("%03d" % t2.julday)
        entry = os.path.join(sdspath, station.net, station.sta, station.comp
                             + '.D', streamData)
        st = read(entry, format="MSEED", starttime=tw['start'],
                  endtime=tw['end'], nearest_sample=True)

        if len(st.get_gaps()) > 0:
            st.merge(method=0, fill_value='interpolate',
                     interpolation_samples=0)

        stream = self.filterWaveform(st)
        return stream

    def readWaveformsPicker_pyrocko(self, station, tw, Origin, ttime):

        obspy_compat.plant()
        cfg = ConfigObj(dict=self.Config)
        if cfg.quantity() == 'displacement':
            try:
                traces = io.load(self.EventPath+'/data/traces_rotated.mseed')
            except:
                traces = io.load(self.EventPath+'/data/traces_restituted.mseed')
        else:
            traces = io.load(self.EventPath+'/data/traces_velocity.mseed')
        for tr in traces:
            tr_name = str(tr.network+'.'+tr.station+'.'+tr.location+'.'
                                    + tr.channel[:3])
            if tr_name == str(station)[:-2] or tr_name == str(station)[:]:
                traces_station = tr
                es = obspy_compat.to_obspy_trace(traces_station)

                st = obspy.Stream()
                st.extend([es])
                stream = ''

                if station.loc == '--':
                    station.loc = ''

                if len(st.get_gaps()) > 0:
                    st.merge(method=0, fill_value='interpolate',
                             interpolation_samples=0)
                stream = self.filterWaveform(st)

                stream.trim(tw['xcorrstart'], tw['xcorrend'])
                return stream

        else:
            pass

    def readWaveformsPicker_colos(self, station, tw, Origin, ttime):

        obspy_compat.plant()
        Config = self.Config
        cfg = ConfigObj(dict=Config)

        traces = io.load(cfg.colosseo_scenario_yml()[:-12]+'scenario.mseed')

        for tr in traces:
            tr_name = str(tr.network+'.'+tr.station+'.'+tr.location+'.'
                          + tr.channel[:3])
        for tr in traces:
            tr_name = str(tr.network+'.'+tr.station+'.'+tr.location+'.'
                                    + tr.channel[:3])
            if tr_name == str(station):
                traces_station = tr

                es = obspy_compat.to_obspy_trace(traces_station)

                st = obspy.Stream()
                st.extend([es])
                stream = ''

                if station.loc == '--':
                    station.loc = ''

                if len(st.get_gaps()) > 0:
                    st.merge(method=0, fill_value='interpolate',
                             interpolation_samples=0)
                stream = self.filterWaveform(st)

                stream.trim(tw['xcorrstart'], tw['xcorrend'])
                return stream

        else:
            pass

    def searchMeta(self, sname, Metalist):

        for i in Metalist:
            if sname == i.getName()[:-2] or sname == i.getName()[:]:
                return i

    def refTrigger(self, RefWaveform, phase):
        Config = self.Config
        cfg = ConfigObj(dict=Config)
        name = ('%s.%s.%s.%s') % (RefWaveform[0].stats.network,
                                  RefWaveform[0].stats.station,
                                  RefWaveform[0].stats.location,
                                  RefWaveform[0].stats.channel)

        i = self.searchMeta(name, self.StationMeta)
        de = loc2degrees(self.Origin, i)

        ptime = 0

        Phase = cake.PhaseDef(phase)
        model = cake.load_model()
        if cfg.colesseo_input() is True:
            arrivals = model.arrivals([de, de], phases=Phase,
                                      zstart=self.Origin.depth, zstop=0.)
        else:
            arrivals = model.arrivals([de, de], phases=Phase,
                                      zstart=self.Origin.depth*km, zstop=0.)
        try:
            ptime = arrivals[0].t
        except Exception:
            arrivals = model.arrivals([de, de], phases=Phase,
                                      zstart=self.Origin.depth*km-0.1)
            ptime = arrivals[0].t

        if ptime == 0:
                raise Exception("\033[31mILLEGAL: phase definition\033[0m")

        tw = self.calculateTimeWindows(ptime)

        if cfg.pyrocko_download() is True:
            stP = self.readWaveformsPicker_pyrocko(i, tw, self.Origin, ptime)
        elif cfg.colesseo_input() is True:
            stP = self.readWaveformsPicker_colos(i, tw, self.Origin, ptime)
        else:
            stP = self.readWaveformsPicker(i, tw, self.Origin, ptime)

        refuntouchname = os.path.basename(self.AF)+'-refstation-raw.mseed'
        stP.write(os.path.join(self.EventPath, refuntouchname), format='MSEED',
                                                                byteorder='>')
        stP.filter("bandpass",
                   freqmin=float(self.Config['refstationfreqmin']),
                   freqmax=float(self.Config['refstationfreqmax']))

        stP.trim(tw['xcorrstart'], tw['xcorrend'])
        trP = stP[0]

        trP.stats.starttime = UTCDateTime(3600)
        refname = os.path.basename(self.AF)+'-refstation-filtered.mseed'
        trP.write(os.path.join(self.EventPath, refname), format='MSEED',
                                                         byteorder='>')

        sta = float(self.Config['refsta'])
        lta = float(self.Config['reflta'])
        cft = recSTALTA(trP.data, int(sta * trP.stats.sampling_rate),
                        int(lta * trP.stats.sampling_rate))

        t = triggerOnset(cft, lta, sta)

        try:
            onset = t[0][0] / trP.stats.sampling_rate

        except Exception:
            onset = self.mintforerun

        trigger = trP.stats.starttime+onset

        tdiff = (trP.stats.starttime + onset)-(UTCDateTime(3600)
                                               + self.mintforerun)

        refp = UTCDateTime(self.Origin.time)+ptime
        reftriggeronset = refp+onset-self.mintforerun

        if int(self.Config['autoxcorrcorrectur']) == 1:
                refmarkername = os.path.join(self.EventPath,
                                             ('%s-marker') % (os.path.basename(
                                              self.AF)))
                fobjrefmarkername = open(refmarkername, 'w')
                fobjrefmarkername.write('# Snuffler Markers File Version\
                                         0.2\n')
                fobjrefmarkername.write(('phase: %s 0 %s    None           None         None         XWStart        None False\n') % (tw['xcorrstart'].strftime('%Y-%m-%d %H:%M:%S.%f'), name))
                fobjrefmarkername.write(('phase: %s 0 %s    None           None         None         XWEnd        None False\n') % (tw['xcorrend'].strftime('%Y-%m-%d %H:%M:%S.%f'), name))
                fobjrefmarkername.write(('phase: %s 1 %s    None           None         None         TheoP        None False\n') % (refp.strftime('%Y-%m-%d %H:%M:%S.%f'), name))
                fobjrefmarkername.write(('phase: %s 3 %s    None           None         None         XTrig        None False') % (reftriggeronset.strftime('%Y-%m-%d %H:%M:%S.%f'), name))
                fobjrefmarkername.close()

                cmd = 'snuffler %s --markers=%s&' % (os.path.join(
                                                    self.EventPath,
                                                    refuntouchname),
                                                    refmarkername)
                os.system(cmd)

                thrOn = float(self.Config['reflta'])
                thrOff = float(self.Config['refsta'])
                plotTrigger(trP, cft, thrOn, thrOff)

                selection = float(input('Enter self picked phase in seconds: '))
                tdiff = selection-self.mintforerun
                refname = os.path.basename(self.AF)+'-shift.mseed'
                trP.stats.starttime = trP.stats.starttime - selection
                trP.write(os.path.join(self.EventPath, refname),
                                       format='MSEED')

        '''
        tdiff = 0
        trigger = trP.stats.starttime
        '''
        To = Trigger(name, trigger, os.path.basename(self.AF), tdiff)

        return tdiff, To

    def shiftSeismograms(self, StreamDict, XcorrDict, pickerShift):

        L = []
        S = []

        dsfactor = float(self.Config['xcorrtreshold'])

        for stream in StreamDict.keys():
            for shift in XcorrDict.keys():
                if stream == shift:

                    StreamDict[stream][0].stats.starttime = StreamDict[stream][0].stats.starttime\
                                                            + XcorrDict[shift].shift+pickerShift

                    if XcorrDict[shift].value > dsfactor:
                        fname = stream + '.new'
                        StreamDict[stream].write(os.path.join(self.AF, fname),
                                                 format='MSEED')

                        if XcorrDict[shift].value < 0:
                            t = -1
                        else:
                            t = 1

                        info = [stream, XcorrDict[shift].shift, t]

                        L.append(info)
                        S.append(stream)
        return L, S

    def writeShift(self, ShiftList):
        import csv
        filename = os.path.join(self.AF, 'shift.dat')
        with open(filename, 'wb') as csv_file:
            writer = csv.writer(csv_file)
            for val in ShiftList:
                writer.writerow([val])

    def f6(self, d1):
        return max(d1, key=d1.get)

    def filterCorrDict(self, CorrDict, onset):

        fCD = {}

        dsfactor = float(self.Config['xcorrtreshold'])
        syn_test = int(self.Config['synthetic_test'])

        if syn_test == 1:
            for stream in CorrDict.keys():
                fCD[stream] = CorrDict[stream]
                fCD[stream].value = fCD[stream].value

        else:
            for stream in CorrDict.keys():
                if abs(CorrDict[stream].value) >= dsfactor:
                    fCD[stream] = CorrDict[stream]
                    fCD[stream].value = fCD[stream].value

        return fCD

    def doXcorr(self, phase, traces):
        StreamDict, SNRDict = self.traveltimes(phase, traces)
        t = self.f6(SNRDict)
        Logfile.add('doXcorr: REFERENCE: ' + t)

        for i in SNRDict.keys():
            Logfile.add('doXcorr: STREAM: ' + i + ' SNR: ' + str(SNRDict[i]))

        alternativeref = os.path.join(*self.AF.split(os.sep)[-1:])+'refstation'

        if self.Config[alternativeref] == '':
            t = t
        else:
            t = self.Config[alternativeref]

        corrDict = {}
        try:
            ref = StreamDict[t][0].data
        except Exception:
            ref = StreamDict[t].data
        Logfile.red('Reference Station of %s for Xcorr Procedure %s'
                    % (os.path.basename(self.AF), t))
        Logfile.red('Enter Xcorr Procedure ')
        for stream in StreamDict.keys():
            a, b = obspy.signal.cross_correlation.xcorr(ref,
                                                        StreamDict[stream][0],
                                                        0)
            shift = a / StreamDict[stream][0].stats.sampling_rate
            corrDict[stream] = Corr(shift, b, a)
            corrDict[stream].value = abs(corrDict[stream].value)
            msg = 'Index: ' + str(a) + ' Value: ' + str(b) + ' ----> '
            msg += (str(stream) + str(StreamDict[stream][0].stats.sampling_rate)
                                + ' SHIFT IN TIME: ' + str(shift))
            Logfile.add(msg)

        Logfile.red('Finish Xcorr Procedure ')
        return corrDict, StreamDict[t], StreamDict

    def doXcorr_dummy(self, phase):
        StreamDict, SNRDict = self.traveltimes(phase)
        t = self.f6(SNRDict)

        alternativeref = os.path.join(*self.AF.split(os.sep)[-1:])+'refstation'

        if self.Config[alternativeref] == '':
            t = t
        else:
            t = self.Config[alternativeref]

        corrDict = {}

        for stream in StreamDict.keys():

            corrDict[stream] = StreamDict[stream]
            corrDict[stream].value = abs(corrDict[stream].value)

        return corrDict, StreamDict[t], StreamDict

    def runXcorr(self, phase, traces):

        CD, ref, WD = self.doXcorr(phase, traces)
        onset = 0
        tdiff, triggerobject = self.refTrigger(ref, phase)

        fCD = self.filterCorrDict(CD, onset)

        return fCD, triggerobject

    def runXcorr_dummy(self, phase):

        CD, ref, WD = self.doXcorr_dummy(phase)
        onset = 0

        fCD = self.filterCorrDict(CD, onset)

        return fCD
