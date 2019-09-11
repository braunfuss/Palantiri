import os
import sys
from pyrocko import guts
import logging
import shutil
import time
import multiprocessing
from optparse import OptionParser
from palantiri.process import optim
from palantiri.common import Basic, Globals, Logfile
from palantiri.common.Program import MainObj
from palantiri.common.ConfigFile import ConfigObj, FilterCfg, OriginCfg, SynthCfg
from collections import OrderedDict
from palantiri.tools import config
from palantiri.tools.config import Event
from palantiri.process import ttt, sembCalc, waveform, times, deserializer
from palantiri.process.array_crosscorrelation_v4 import Xcorr, cmpFilterMetavsXCORR
from pyrocko import util, io

import numpy as num
if sys.version_info.major >= 3:
    import _pickle as pickle
    xrange = range
else:
    import cPickle as pickle


rstate = num.random.RandomState()
logger = logging.getLogger(sys.argv[0])
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(message)s")
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)
evpath = None


def initModule():

    global evpath
    global evpath_emp

    parser = OptionParser(usage="%prog -f Eventpath -e Eventpath_smaller(optional)")
    parser.add_option("-f", "--evpath", type="string", dest="evpath",
                      help="evpath")
    parser.add_option("-e", "--evpath_emp", type="string", dest="evpath_emp",
                      help="evpath_emp")

    (options, args) = parser.parse_args()
    if options.evpath is None:
        parser.error("non existing eventpath")
        return False

    evpath = options.evpath
    evpath_emp = options.evpath_emp
    Globals.setEventDir(evpath)
    Globals.setEventDir_emp(evpath_emp)
    return True


def processLoop(traces=None, stations=None, cluster=None):

    C = config.Config(evpath, eventpath_emp=evpath_emp)
    Origin = C.parseConfig('origin')
    flag_rpe = False

    try:
        Syn_in = C.parseConfig('syn')
        syn_in = SynthCfg(Syn_in)
    except TypeError:
        pass
    try:
        Syn_in_emp = C.parseConfig('syn_emp')
        syn_in_emp = SynthCfg(Syn_in_emp)
    except IndexError:
        pass
    Config = C.parseConfig('config')

    cfg = ConfigObj(dict=Config)
    phases = cfg.Str('ttphases')
    phases = phases.split(',')

    if cfg.pyrocko_download() is True:
        Meta = C.readpyrockostations(phases)
    elif cfg.colesseo_input() is True:
        scenario = guts.load(filename=cfg.colosseo_scenario_yml())
        scenario_path = cfg.colosseo_scenario_yml()[:-12]
        Meta = C.readcolosseostations(scenario_path)
    else:
        Meta = C.readMetaInfoFile()
    Folder = C.createFolder()
    C.writeConfig(Config, Origin, Folder)

    filter = FilterCfg(Config)
    if cfg.UInt('forerun') > 0:
        ntimes = int((cfg.UInt('forerun') + cfg.UInt('duration')) /
                     cfg.Float('step'))
    else:
        ntimes = int((cfg.UInt('duration')) / cfg.Float('step'))
    origin = OriginCfg(Origin)

    if cfg.colesseo_input() is True:
        events = scenario.get_events()
        ev = events[0]
        origin.strike = str(ev.moment_tensor.strike1)
        origin.rake = str(ev.moment_tensor.rake1)
        origin.dip = str(ev.moment_tensor.dip1)
        strike = ev.moment_tensor.strike1
        origin.lat = str(ev.lat)
        origin.lon = str(ev.lon)
        origin.depth = str(ev.depth/1000.)
        depth = ev.depth
        origin.time = util.time_to_str(ev.time)
        time_ev = util.time_to_str(ev.time)
        lat = ev.lat
        lon = ev.lon
        rake = ev.moment_tensor.rake1
        dip = ev.moment_tensor.dip1
        Origin['strike'] = str(ev.moment_tensor.strike1)
        Origin['rake'] = str(ev.moment_tensor.rake1)
        Origin['dip'] = str(ev.moment_tensor.dip1)
        Origin['lat'] = str(ev.lat)
        Origin['lon'] = str(ev.lon)
        Origin['time'] = util.time_to_str(ev.time)
        Origin['depth'] = str(ev.depth/1000.)
        ev = Event(lat, lon, depth, time_ev,
                   strike=strike, dip=dip, rake=rake)
    else:

        default = 0
        strike = origin.strike(default)
        dip = origin.dip(default)
        rake = origin.rake(default)

        ev = Event(origin.lat(), origin.lon(), origin.depth(), origin.time(),
                   strike=strike, dip=dip, rake=rake)

    if cfg.Bool('correct_shifts_empirical') is True:

        Origin_emp = C.parseConfig('origin_emp')
        origin_emp = OriginCfg(Origin_emp)
        ev_emp = Event(origin_emp.lat(), origin_emp.lon(), origin_emp.depth(),
                       origin_emp.time(), strike=strike, dip=dip, rake=rake)
    filtername = filter.filterName()
    Logfile.add('filtername = ' + filtername)

    XDict = OrderedDict()
    RefDict = OrderedDict()
    SL = OrderedDict()
    refshifts_global = []
    newFreq = str(filter.newFrequency())
    xcorrnetworks = cfg.String('networks').split(',')

    if cfg.Int('xcorr') is 1:

        fobjreferenceshiftname = newFreq + '_' + filtername + '.refpkl'
        rp = os.path.join(Folder['semb'], fobjreferenceshiftname)
        fobjreferenceshiftnameemp = newFreq + '_' + filtername + 'emp' + '.refpkl'
        rpe = os.path.join(Folder['semb'], fobjreferenceshiftnameemp)
        fobjpickleshiftname = newFreq + '_' + filtername + '.xcorrpkl'
        ps = os.path.join(Folder['semb'], fobjpickleshiftname)
        if cfg.Bool('synthetic_test') is False:
            if cfg.quantity() == 'displacement':
                try:
                    traces = io.load(evpath+'/data/traces_rotated.mseed')
                except Exception:
                    traces = io.load(evpath+'/data/traces_restituted.mseed')
            else:
                traces = io.load(evpath+'/data/traces_velocity.mseed')
        if(os.path.isfile(rp) and os.path.getsize(rp) != 0
           and os.path.isfile(ps) and os.path.getsize(ps) != 0):
            Logfile.add('xcorr/reference shift file exits : ' + rp)
            Logfile.add('loaded reference shift')

            if sys.version_info.major >= 3:
                f = open(rp, 'rb')
            else:
                f = open(rp)

            RefDict = pickle.load(f)
            if sys.version_info.major >= 3:
                x = open(ps, 'rb')
            else:
                x = open(ps)
            XDict = pickle.load(x)
            for i in xcorrnetworks:
                SL[i] = len(Config[i].split('|'))
        else:
            if cfg.Bool('synthetic_test') is False:
                if cfg.quantity() == 'displacement':
                    try:
                        traces = io.load(evpath+'/data/traces_rotated.mseed')
                    except Exception:
                        traces = io.load(evpath+'/data/traces_restituted.mseed')
                else:
                    traces = io.load(evpath+'/data/traces_velocity.mseed')
            SL = {}
            for i in xcorrnetworks:
                W = {}
                network = cfg.String(i).split('|')
                FilterMeta = ttt.filterStations(Meta, Config, Origin, network)
                arrayfolder = os.path.join(Folder['semb'], i)

                if os.access(arrayfolder, os.F_OK) is False:
                    os.makedirs(arrayfolder)
                if cfg.pyrocko_download() is True:
                    # TODO check seperate xcoor nescessity
                    A = Xcorr(ev, FilterMeta, evpath, Config, Syn_in,
                              arrayfolder)
                print("run Xcorr")
                phase = phases[0]
                W, triggerobject = A.runXcorr(phase, traces)
                XDict[i] = W
                RefDict[i] = triggerobject.tdiff
                SL[i] = len(network)
                for j in range(0, len(FilterMeta)):
                    refshifts_global.append(triggerobject.tdiff)

            if sys.version_info.major >= 3:
                fobjrefshift = open(rp, 'wb')
            else:
                fobjrefshift = open(rp, 'w')
            pickle.dump(RefDict, fobjrefshift)
            fobjrefshift.close()

            if sys.version_info.major >= 3:
                output = open(ps, 'wb')
            else:
                output = open(ps, 'w')
            pickle.dump(XDict, output)
            output.close()

    else:
        fobjreferenceshiftname = newFreq + '_' + filtername + '.refpkl'
        rp = os.path.join(Folder['semb'], fobjreferenceshiftname)
        fobjreferenceshiftnameemp = newFreq + '_' + filtername + 'emp' + '.refpkl'
        rpe = os.path.join(Folder['semb'], fobjreferenceshiftnameemp)
        fobjpickleshiftname = newFreq + '_' + filtername + '.xcorrpkl'
        ps = os.path.join(Folder['semb'], fobjpickleshiftname)
        refshift = 0
        if(os.path.isfile(rp) and os.path.getsize(rp) != 0
           and os.path.isfile(ps) and os.path.getsize(ps) != 0):
            Logfile.add('Temporay Memory file exits : ' + rp)
            if sys.version_info.major >= 3:
                f = open(rp, 'rb')
            else:
                f = open(rp)

            RefDict = pickle.load(f)
            if sys.version_info.major >= 3:
                x = open(ps, 'rb')
            else:
                x = open(ps)
            XDict = pickle.load(x)

            for i in xcorrnetworks:
                SL[i] = len(Config[j].split('|'))
                network = cfg.String(i).split('|')
                FilterMeta = ttt.filterStations(Meta, Config, Origin, network)
                RefDict[i] = refshift

                for j in range(0, len(FilterMeta)):
                    refshifts_global.append(refshift)
        else:
            SL = {}
            for i in xcorrnetworks:
                W = {}
                refshift = 0
                network = cfg.String(i).split('|')
                FilterMeta = ttt.filterStations(Meta, Config, Origin, network)
                arrayfolder = os.path.join(Folder['semb'], i)

                if os.access(arrayfolder, os.F_OK) is False:
                    os.makedirs(arrayfolder)
                if cfg.pyrocko_download() is True:
                    # TODO check seperate xcoor nescessity
                    A = Xcorr(ev, FilterMeta, evpath, Config, Syn_in,
                              arrayfolder)
                else:
                    A = Xcorr(ev, FilterMeta, evpath, Config, Syn_in,
                              arrayfolder)

                print("run Xcorr")
                phase = phases[0]
                W, triggerobject = A.runXcorr_dummy(phase)

                XDict[j] = W
                RefDict[j] = refshift
                SL[j] = len(network)
                for j in range(0, len(FilterMeta)):
                    refshifts_global.append(refshift)

            if sys.version_info.major >= 3:
                fobjrefshift = open(rp, 'wb')
            else:
                fobjrefshift = open(rp, 'w')
            pickle.dump(RefDict, fobjrefshift)
            fobjrefshift.close()

            if sys.version_info.major >= 3:
                output = open(ps, 'wb')
            else:
                output = open(ps, 'w')
            pickle.dump(XDict, output)
            output.close()

    if sys.version_info.major >= 3:
        for j in sorted(XDict.keys()):
            Logfile.red('Array %s has %3d of %3d Stations left' %
                        (j, len(XDict[j]), SL[j]))
    else:
        for j in sorted(XDict.keys()):
            Logfile.red('Array %s has %3d of %3d Stations left' %
                        (j, len(XDict[j]), SL[j]))
    while True:
        if sys.version_info.major >= 3:
            nnl = input("please enter your choice: ")
        else:
            nnl = raw_input("please enter your choice: ")

        if len(nnl) == 0:
            if not Basic.question('Process all networks ?'):
                continue

            Logfile.red('This networks will be used for processing: %s'
                        % (Config['networks']))
            break

        elif str(nnl) == 'quit':
            sys.exit()

        elif str(nnl) == 'rerun':
            event = os.path.join(*evpath.split('/')[-1:])

            try:
                os.remove(rp)
                os.remove(ps)

            except Exception:
                pass

            mainfolder = os.path.join(os.path.sep, *evpath.split('/')[:-2])
            os.chdir(mainfolder)

            cmd = ('%s arraytool.py process %s') % (sys.executable, event)
            Logfile.add('cmd = ' + cmd)
            os.system(cmd)
            sys.exit()

        else:

            names = nnl.split(',')
            isOk = True

            for array in names:
                arrayfolder = os.path.join(Folder['semb'], array)

                if not os.path.isdir(arrayfolder):
                    Logfile.error('Illegal network name ' + str(array))
                    isOk = False
                    break
            if not isOk:
                continue

            Logfile.add('This networks will be used for processing: %s'
                        % (nnl))
            Config['networks'] = nnl
            break

    for j in range(3, 0, -1):
        time.sleep(1)
        Logfile.red('Start processing in %d seconds ' % (j))

    wd = Origin['depth']
    start, stop, step = cfg.String('depths').split(',')

    start = int(start)
    stop = int(stop)+1
    step_depth = int(step)
    filters = cfg.String('filters')
    filters = int(filters)
    Logfile.add('working on ' + Config['networks'])
    if cfg.Bool('correct_shifts_empirical') is True:
        emp_loop = True
    else:
        emp_loop = False

# ==================================loop over phases======================
    for phase in phases:
        if phase is 'P':
            desired = 'Z'
        if phase is 'S':
            desired = 'T'
        # ==================================loop over filter setups=====
        for filterindex in xrange(0, filters):
            # ==================================loop over depth=======
            for depthindex in xrange(start, stop, step_depth):

                workdepth = float(wd) + depthindex
                if cfg.Int('dimz') == 0:
                    Origin['depth'] = workdepth
                ev = Event(Origin['lat'], Origin['lon'], Origin['depth'],
                           Origin['time'], strike=strike, dip=dip, rake=rake)
                Logfile.add('WORKDEPTH: ' + str(Origin['depth']))
                networks = Config['networks'].split(',')

                ASL = []
                weights = []
                array_centers = []
                counter = 1
                stations_per_array = []
                Wdfs = []
                Wdfs_emp = []
                FilterMetas = []
                TTTgrids = OrderedDict()
                mints = []
                maxts = []
                for i in networks:
                    refshifts = []

                    arrayname = i
                    arrayfolder = os.path.join(Folder['semb'], arrayname)

                    network = Config[i].split('|')
                    Logfile.add('network: ' + str(network))

                    FilterMeta = ttt.filterStations(Meta, Config, Origin,
                                                    network)

                    W = XDict[i]
                    refshift = RefDict[i]
                    for j in range(0, len(FilterMeta)):
                        if cfg.correct_shifts() is False:
                            refshift = refshift*0.
                        refshifts.append(refshift)

                    FilterMeta = cmpFilterMetavsXCORR(W, FilterMeta)

                    Logfile.add('BOUNDING BOX DIMX: %s  DIMY: %s  GRIDSPACING:\
                                %s \n' % (Config['dimx'], Config['dimy'],
                                          Config['gridspacing']))

                    Logfile.red('Calculating Traveltime Grid')
                    t1 = time.time()

                    isParallel = False
                    TTTGridMap = []
                    mint = []
                    maxt = []
                    ttt_model = cfg.Str('traveltime_model')
                    try:
                        px = os.path.abspath(os.path.join(os.getcwd(),
                                             os.pardir))
                        pdx = str('/tttgrid/tttgrid%s_%s_%s_%s_%s.pkl'
                                  % (phase, ttt_model, ev.time, arrayname,
                                     workdepth))
                        px_path = px+pdx
                        f = open(px_path, 'rb')

                        print("loading travel time tttgrid%s_%s_%s_%s_%s.pkl"
                              % (phase, ttt_model, ev.time, arrayname,
                                 workdepth))
                        TTTGridMap, mint, maxt = pickle.load(f)
                        f.close()

                        print("loading of travel time grid sucessful")
                    except Exception:
                        print("loading of travel time grid unsucessful,\n \
                              will now calculate the grid:")
                        if isParallel:
                            maxp = 6
                            po = multiprocessing.Pool(maxp)

                            for i in xrange(len(FilterMeta)):
                                po.apply_async(ttt.calcTTTAdv,
                                               (Config, FilterMeta[i], Origin,
                                                i, arrayname, W, refshift))

                                po.close()
                                po.join()
                        else:
                            for i in xrange(len(FilterMeta)):
                                if cfg.Int('dimz') != 0:
                                    t1 = time.time()
                                    ttt.calcTTTAdv_cube(Config, FilterMeta[i],
                                                        Origin, i, arrayname,
                                                        W, refshift, phase)

                                    Logfile.add('ttt.calcTTTAdv : '
                                                + str(time.time() - t1)
                                                + ' sec.')
                                else:
                                    t1 = time.time()
                                    ttt.calcTTTAdv(Config, FilterMeta[i],
                                                   Origin,
                                                   i, arrayname, W, refshift,
                                                   phase)

                                    Logfile.add('ttt.calcTTTAdv : '
                                                + str(time.time() - t1)
                                                + ' sec.')
                        assert len(FilterMeta) > 0
                        if cfg.Int('dimz') != 0:
                            TTTGridMap = deserializer.deserializeTTT_cube(len(FilterMeta))
                            mint, maxt = deserializer.deserializeMinTMaxT(len(FilterMeta))
                        else:
                            TTTGridMap = deserializer.deserializeTTT(len(FilterMeta))
                            mint, maxt = deserializer.deserializeMinTMaxT(len(FilterMeta))

                        px = os.path.abspath(os.path.join(os.getcwd(),
                                             os.pardir))
                        pdx = str('/tttgrid/tttgrid%s_%s_%s_%s_%s.pkl'
                                 % (phase, ttt_model, ev.time, arrayname,
                                    workdepth))
                        px_path = px+pdx
                        f = open(px_path, 'wb')
                        pickle.dump([TTTGridMap, mint, maxt], f)
                        f.close()
                    t2 = time.time()
                    Logfile.red('%s took %0.3f s' % ('TTT', (t2-t1)))

                    if cfg.Bool('correct_shifts_empirical') is True:

                        Logfile.add('BOUNDING BOX DIMX: %s  DIMY: %s  GRIDSPACING:\
                                    %s \n' % (Config['dimx_emp'],
                                              Config['dimy_emp'],
                                              Config['gridspacing']))

                        try:
                            f = open(os.path.abspath(os.path.join(os.getcwd(),
                                                     os.pardir))+'/tttgrid/tttgrid%s_%s_%s_%s_%s_emp.pkl'
                                     % (phase, ttt_model, ev_emp.time,
                                        arrayname, workdepth), 'rb')
                            print("loading travel time grid%s_%s_%s_%s_%s_emp.pkl"
                                  % (phase, ttt_model, ev_emp.time, arrayname,
                                     workdepth))
                            TTTGridMap_emp, mint_emp, maxt_emp = pickle.load(f)
                            f.close()

                            print("loading of travel time grid sucessful")
                        except Exception:
                            print("loading of travel time grid unsucessful,\n \
                                  will now calculate the grid:")
                            if isParallel:
                                maxp = 6
                                po = multiprocessing.Pool(maxp)

                                for i in xrange(len(FilterMeta)):
                                    po.apply_async(ttt.calcTTTAdv,
                                                   (Config, FilterMeta[i],
                                                    Origin, i, arrayname,
                                                    W, refshift))

                                    po.close()
                                    po.join()
                            else:
                                for i in xrange(len(FilterMeta)):
                                    t1 = time.time()
                                    ttt.calcTTTAdv(Config, FilterMeta[i],
                                                   Origin_emp,
                                                   i, arrayname, W, refshift,
                                                   phase, flag_rpe=True)

                                    assert len(FilterMeta) > 0
                                    TTTGridMap_emp = deserializer.deserializeTTT(len(FilterMeta))
                                    mint_emp, maxt_emp = deserializer.deserializeMinTMaxT(len(FilterMeta))
                                    f = open(os.path.abspath(os.path.join(os.getcwd(),
                                             os.pardir))+'/tttgrid/tttgrid%s_\
                                             %s_%s_%s_%s_emp.pkl'
                                             % (phase, ttt_model, ev_emp.time,
                                                arrayname, workdepth), 'wb')
                                    print("dumping the traveltime grid for this array")
                                    pickle.dump([TTTGridMap_emp, mint_emp,
                                                 maxt_emp], f)
                                    f.close()
                        t2 = time.time()
                        Logfile.red('%s took %0.3f s' % ('TTT', (t2-t1)))

                    switch = filterindex
                    tw = times.calculateTimeWindows(mint, maxt, Config,
                                                    ev, switch)
                    if cfg.Bool('correct_shifts_empirical') is True:
                        tw_emp = times.calculateTimeWindows(mint_emp, maxt_emp,
                                                            Config, ev_emp,
                                                            switch)

                        if cfg.pyrocko_download() is True:
                            if cfg.quantity() == 'displacement':
                                Wd_emp = waveform.readWaveformsPyrocko_restituted(
                                    FilterMeta, tw, evpath, ev_emp, desired)
                            elif cfg.Bool('synthetic_test') is True:
                                Wd_emp = waveform.readWaveformsPyrockodummy(FilterMeta,
                                                                        tw_emp,
                                                                        evpath_emp,
                                                                        ev_emp)
                            else:
                                Wd_emp = waveform.readWaveformsPyrocko(FilterMeta,
                                                                       tw_emp,
                                                                       evpath_emp,
                                                                       ev_emp,
                                                                       desired)
                        elif cfg.colesseo_input() is True:
                            Wd_emp = waveform.readWaveforms_colesseo(FilterMeta,
                                                                     tw_emp,
                                                                     evpath_emp,
                                                                     ev_emp, C)
                        else:
                            Wd_emp = waveform.readWaveforms(FilterMeta, tw_emp,
                                                            evpath_emp, ev_emp)
                        if cfg.Bool('synthetic_test') is True\
                           or cfg.Bool('dynamic_filter') is True:
                            Wdf_emp = waveform.processdummyWaveforms(Wd_emp, Config,
                                                                 Folder, arrayname,
                                                                 FilterMeta, ev_emp,
                                                                 switch, W)
                            Wdfs_emp.extend(Wdf_emp)
                        else:
                            if switch == 0:
                                ff1 = filter.flo()
                                ff2 = filter.fhi()
                            if switch == 1:
                                ff1 = filter.flo2()
                                ff2 = filter.fhi2()
                            ps_wdf_emp = os.path.join(Folder['semb'],
                                                      "fobjpickle_process_emp\
                                                      _%s_%s%s"
                                                      % (arrayname, ff1, ff2))
                            if cfg.Bool('load_wdf') is True:
                                try:
                                    f = open(ps_wdf_emp, 'rb')
                                    Wdf_emp = pickle.load(f)
                                except Exception:
                                    Wdf_emp = waveform.processWaveforms(Wd_emp,
                                                                        Config,
                                                                        Folder,
                                                                        arrayname,
                                                                        FilterMeta,
                                                                        ev_emp,
                                                                        switch,
                                                                        W)

                                    fobj_proc = open(ps_wdf_emp, 'wb')
                                    pickle.dump(Wdf_emp, fobj_proc)
                                    f = open(ps_wdf_emp, 'rb')
                                    Wdf_emp = pickle.load(f)
                            else:
                                Wdf_emp = waveform.processWaveforms(Wd_emp,
                                                                    Config,
                                                                    Folder,
                                                                    arrayname,
                                                                    FilterMeta,
                                                                    ev_emp,
                                                                    switch, W)
                            Wdfs_emp.extend(Wdf_emp)
                    if cfg.pyrocko_download() is True:
                        if cfg.quantity() == 'displacement':
                            Wd = waveform.readWaveformsPyrocko_restituted(
                                FilterMeta, tw, evpath, ev, desired)
                        elif cfg.Bool('synthetic_test') is True:
                            Wd = waveform.readWaveformsPyrockodummy(FilterMeta,
                                                                    tw, evpath,
                                                                    ev)
                        else:
                            Wd = waveform.readWaveformsPyrocko(FilterMeta, tw,
                                                               evpath, ev,
                                                               desired)
                    elif cfg.colesseo_input() is True:
                        Wd = waveform.readWaveforms_colesseo(FilterMeta, tw,
                                                             evpath, ev, C)
                    else:
                        Wd = waveform.readWaveforms(FilterMeta, tw, evpath, ev)
                    if cfg.Bool('synthetic_test') is True\
                       or cfg.Bool('dynamic_filter') is True:
                        Wdf = waveform.processdummyWaveforms(Wd, Config,
                                                             Folder, arrayname,
                                                             FilterMeta, ev,
                                                             switch, W)
                        Wdfs.extend(Wdf)
                    else:
                        if switch == 0:
                            ff1 = filter.flo()
                            ff2 = filter.fhi()
                        if switch == 1:
                            ff1 = filter.flo2()
                            ff2 = filter.fhi2()
                        ps_wdf = os.path.join(Folder['semb'], "fobjpickle_process_%s_%s%s" % (arrayname, ff1, ff2))
                        if cfg.Bool('load_wdf') is True:
                            try:
                                f = open(ps_wdf, 'rb')
                                Wdf = pickle.load(f)
                            except:
                                Wdf = waveform.processWaveforms(Wd, Config,
                                                                Folder,
                                                                arrayname,
                                                                FilterMeta,
                                                                ev, switch, W)

                                fobj_proc = open(ps_wdf, 'wb')
                                pickle.dump(Wdf, fobj_proc)
                                f = open(ps_wdf, 'rb')
                                Wdf = pickle.load(f)
                        else:
                            Wdf = waveform.processWaveforms(Wd, Config, Folder,
                                                            arrayname,
                                                            FilterMeta,
                                                            ev, switch, W)
                        Wdfs.extend(Wdf)

                    C.writeStationFile(FilterMeta, Folder, counter)
                    Logfile.red('%d Streams added for Processing' % (len(Wd)))

                    t1 = time.time()

                    f = open(os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/tttgrid/tttgrid%s_%s_%s_%s_%s.pkl'
                             % (phase, ttt_model, ev.time, arrayname,
                                workdepth), 'rb')
                    TTTGridMap, mint, maxt = pickle.load(f)
                    f.close()
                    if switch == 0:
                        step = cfg.step()
                    if switch == 1:
                        step = cfg.step_f2()
                    if cfg.UInt('forerun') > 0:
                        ntimes = int((cfg.UInt('forerun') +
                                      cfg.UInt('duration')) / step)
                    else:
                        ntimes = int((cfg.UInt('duration')) / step)
                    if cfg.Bool('combine_all') is False:

                        if cfg.optimize() is True:
                            optim.solve(counter, Config, Wdf, FilterMeta,
                                        mint, maxt, TTTGridMap, Folder,
                                        Origin, ntimes, switch, ev,
                                        arrayfolder, syn_in, refshifts, phase,
                                        rpe+str(arrayname), flag_rpe)
                        else:
                            if cfg.Bool('correct_shifts_empirical') is True:
                                if cfg.Bool('correct_shifts_empirical_run') is True:
                                    winlen_emp = cfg.winlen_emp()
                                    step_emp = cfg.step_emp()
                                    if cfg.UInt('forerun') > 0:
                                        ntimes_emp = int((cfg.UInt('forerun_emp') + cfg.UInt('duration_emp'))/step_emp)
                                    else:
                                        ntimes_emp = int((cfg.UInt('duration_emp')) / step_emp)
                                        f = open(os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/tttgrid/tttgrid%s_%s_%s_%s_%s_emp.pkl'
                                             % (phase, ttt_model, ev_emp.time,
                                                arrayname,
                                                workdepth), 'rb')
                                    TTTGridMap_emp, mint_emp, maxt_emp = pickle.load(f)
                                    f.close()
                                    flag_rpe = True
                                    arraySemb, weight, array_center = sembCalc.doCalc(
                                        counter, Config, Wdf_emp, FilterMeta,
                                        mint_emp, maxt_emp,
                                        TTTGridMap_emp, Folder, Origin_emp,
                                        ntimes_emp, switch, ev_emp,
                                        arrayfolder, syn_in_emp, refshifts,
                                        phase, rpe+str(arrayname), flag_rpe)
                                    if sys.version_info.major >= 3:
                                        f = open(rpe+str(arrayname), 'rb')
                                    else:
                                        f = open(rpe+str(arrayname))
                                    RefDict_empirical = pickle.load(f)
                                    refshifts = RefDict_empirical
                                    for j in range(0, len(FilterMeta)):
                                        if cfg.correct_shifts() is False:
                                            refshifts[j] = refshifts[j]*0.
                            flag_rpe = False
                            f = open(os.path.abspath(os.path.join(os.getcwd(),
                                     os.pardir))+'/tttgrid/tttgrid%s_%s_%s_%s_%s.pkl'
                                     % (phase, ttt_model, ev.time, arrayname,
                                        workdepth), 'rb')
                            print("loading travel time grid%s_%s_%s_%s_%s.pkl"
                                  % (phase, ttt_model, ev.time, arrayname,
                                      workdepth))
                            TTTGridMap_mew, mintt, maxtt = pickle.load(f)
                            f.close()

                            arraySemb, weight, array_center = sembCalc.doCalc(
                                counter, Config, Wdf, FilterMeta, mintt, maxtt,
                                TTTGridMap_mew, Folder, Origin, ntimes, switch,
                                ev, arrayfolder, syn_in, refshifts, phase,
                                rpe+str(arrayname), flag_rpe, len(FilterMeta))
                            weights.append(weight)
                            array_centers.append(array_center)
                            ASL.append(arraySemb)
                            sembCalc.writeSembMatricesSingleArray(arraySemb,
                                                                  Config,
                                                                  Origin,
                                                                  arrayfolder,
                                                                  ntimes,
                                                                  switch,
                                                                  phase)

                    fileName = os.path.join(arrayfolder, 'stations.txt')
                    Logfile.add('Write to file ' + fileName)

                    fobjarraynetwork = open(fileName, 'w')
                    for i in FilterMeta:
                        fobjarraynetwork.write(('%s %s %s\n') %
                                               (i.getName(), i.lat, i.lon))

                    fobjarraynetwork.close()
                    t2 = time.time()
                    Logfile.add('CALC took %0.3f sec' % (t2-t1))
                    counter += 1
                    stations_per_array.append(len(FilterMeta))
                    TTTgrids.update(TTTGridMap)
                    mints.append(mint)
                    maxts.append(maxt)
                    FilterMetas[len(FilterMetas):] = FilterMeta
                    TTTGridMap = []

                if cfg.Bool('combine_all') is True:
                    if cfg.pyrocko_download() is True:
                        if cfg.Bool('synthetic_test') is True:
                            Wd = waveform.readWaveformsPyrockodummy(
                                    FilterMetas, tw, evpath, ev)
                        else:
                            if cfg.quantity() == 'displacement':
                                Wd = waveform.readWaveformsPyrocko_restituted(
                                    FilterMetas, tw, evpath, ev, desired)
                            else:
                                Wd = waveform.readWaveformsPyrocko(FilterMetas,
                                                                   tw, evpath,
                                                                   ev, desired)
                    elif cfg.colesseo_input() is True:
                        Wd = waveform.readWaveforms_colesseo(FilterMetas, tw,
                                                             evpath, ev, C)
                    else:
                        Wd = waveform.readWaveforms(FilterMetas, tw, evpath,
                                                    ev)
                    if cfg.Bool('synthetic_test') is True:
                        Wdf = waveform.processdummyWaveforms(Wd, Config,
                                                             Folder, arrayname,
                                                             FilterMetas, ev,
                                                             switch, W)
                    else:
                        if switch == 0:
                            ff1 = filter.flo()
                            ff2 = filter.fhi()
                        if switch == 1:
                            ff1 = filter.flo2()
                            ff2 = filter.fhi2()
                        ps_wdf = os.path.join(Folder['semb'],
                                              "fobjpickle_process_%s_%s%s_\
                                              combined" % (arrayname, ff1, ff2))
                        if cfg.Bool('load_wdf') is True:
                            try:
                                f = open(ps_wdf, 'rb')
                                Wdf = pickle.load(f)
                                print('loaded wdf')
                                print(ps_wdf)
                            except:
                                Wdf = waveform.processWaveforms(Wd, Config,
                                                                Folder,
                                                                arrayname,
                                                                FilterMetas,
                                                                ev, switch, W)

                                fobj_proc = open(ps_wdf, 'wb')
                                pickle.dump(Wdf, fobj_proc)
                                print('dumped wdf')
                        else:
                            Wdf = waveform.processWaveforms(Wd, Config, Folder,
                                                            arrayname,
                                                            FilterMetas,
                                                            ev, switch, W)

                    mint = num.min(mints)
                    maxt = num.max(maxts)
                    nstats = stations_per_array
                    flag_rpe = False
                    if cfg.Bool('bootstrap_array_weights') is False:
                        arraySemb, weight, array_center = sembCalc.doCalc(
                            counter, Config, Wdf, FilterMetas, mint, maxt,
                            TTTgrids, Folder, Origin, ntimes, switch,
                            ev, arrayfolder, syn_in, refshifts_global, phase,
                            rpe+str(arrayname), flag_rpe, nstats)
                        ASL.append(arraySemb)
                        weights.append(weight)
                        array_centers.append(array_center)
                        sembCalc.writeSembMatricesSingleArray(arraySemb,
                                                              Config, Origin,
                                                              arrayfolder,
                                                              ntimes, switch,
                                                              phase)
                    else:
                        nboot = cfg.Int('n_bootstrap')
                        tmp_general = 1
                        for ibootstrap in range(nboot):
                            f = rstate.uniform(0., 1., size=counter+1)
                            f = num.sort(f)
                            g = f[1:] - f[:-1]
                            k = 0
                            ws = []

                            for wss in range(0, counter-1):
                                for stats in range(0, stations_per_array[k]):
                                    ws.append(g[k])
                                k =+ 1
                            ws = num.asarray(ws)
                            arraySemb, weight, array_center = sembCalc.doCalc(
                                counter, Config, Wdf, FilterMetas, mint, maxt,
                                TTTgrids, Folder, Origin, ntimes, switch,
                                ev, arrayfolder, syn_in, refshifts_global,
                                phase, rpe+str(arrayname), flag_rpe, nstats,
                                bs_weights=ws)

                            ASL.append(arraySemb)
                            weights.append(weight)
                            array_centers.append(array_center)
                            sembCalc.writeSembMatricesSingleArray(arraySemb,
                                                                  Config,
                                                                  Origin,
                                                                  arrayfolder,
                                                                  ntimes,
                                                                  switch,
                                                                  phase,
                                                                  bootstrap=ibootstrap)

                            if ASL:
                                Logfile.red('collect semblance matrices from\
                                            all arrays')
                                sembmax, tmp = sembCalc.collectSemb(ASL,
                                                                    Config,
                                                                    Origin,
                                                                    Folder,
                                                                    ntimes,
                                                                    len(networks),
                                                                    switch,
                                                                    array_centers,
                                                                    phase,
                                                                    cboot=ibootstrap)
                                tmp_general *= tmp
                                ASL = []
                        sembmax, tmp = sembCalc.collectSemb(ASL, Config,
                                                            Origin, Folder,
                                                            ntimes,
                                                            len(networks),
                                                            switch,
                                                            array_centers,
                                                            phase,
                                                            cboot=None,
                                                            temp_comb=tmp_general)
                if cfg.optimize_all() is True:
                    import optim_csemb
                    sembmax, tmp = sembCalc.collectSemb(ASL, Config, Origin,
                                                        Folder, ntimes,
                                                        len(networks), switch)
                    optim_csemb.solve(counter, Config, Wdf, FilterMeta, mint,
                                      maxt, TTTGridMap, Folder, Origin,
                                      ntimes, switch, ev, arrayfolder,
                                      syn_in, ASL, sembmax, evpath, XDict,
                                      RefDict, workdepth, filterindex, Wdfs)

                if ASL:
                    Logfile.red('collect semblance matrices from all arrays')
                    sembmax, tmp = sembCalc.collectSemb(ASL, Config, Origin,
                                                        Folder,
                                                        ntimes, len(networks),
                                                        switch, array_centers,
                                                        phase)
                    if cfg.Bool('weight_by_noise') is True:
                        sembCalc.collectSembweighted(ASL, Config, Origin,
                                                     Folder, ntimes,
                                                     len(networks), switch,
                                                     weights)

    else:
        Logfile.red('Nothing to do  -> Finish')
    print("last work depth:")
    print(workdepth)


class ProcessMain(MainObj):

    def __init__(self):
        initModule()

        MainObj.__init__(self, self, '0.3', 'process_run.log', 'process.log')

    def init(self):
        return True

    def process(self):
        processLoop()
        return True

    def finish(self):
        pass


def MainProc():

    mainObj = ProcessMain()

    mainObj.run()


isClient = False


def main():

    if not isClient:
        MainProc()
