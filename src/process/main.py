import os
import sys
from pyrocko import guts
import logging
import shutil
import time
import glob
import multiprocessing
from optparse import OptionParser
from palantiri.process import optim
from palantiri.common import Basic, Globals, Logfile
from palantiri.common.Program import MainObj
from palantiri.common.ConfigFile import ConfigObj, FilterCfg, OriginCfg, SynthCfg
from collections import OrderedDict
from palantiri.tools import config

from palantiri.tools.config import Event, Config, ClusterConfig, PalantiriConfig, PalantiriDataConfig, PalantiriXcorrConfig, PalantiriFilterConfig, PalantiriWeightConfig, PalantiriGeometryConfig, PalantiriSyntheticConfig
from palantiri.process import ttt, sembCalc, waveform, times, deserializer
from palantiri.process.array_crosscorrelation_v4 import Xcorr, cmpFilterMetavsXCORR
from pyrocko import util, io, guts
from pyrocko.guts import Object, Float, Int, String, Bool, List, Tuple

import subprocess
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


def check_is_empty(evpath, move=False):
    from pathlib import Path
    folder_used = False
    path = Path(evpath+'/work/semblance/').glob('*.ASC')
    for p in path:
        folder_used = True
    if folder_used is True:
        if move is False:
            print('Workdir is not empty, to overwrite use --force')
            quit()
        else:
            for kiter in range(0, 100):
                new_folder_used = False
                print('copying backup of work -might take some time')
                path = Path(evpath+'/work_%s/semblance/' %kiter).glob('*.ASC')
                for p in path:
                    new_folder_used = True
                if new_folder_used is False:
                    os.system('cp -r %s %s_%s' % (evpath+'/work/semblance', evpath+'/work/semblance', kiter))
                    print(evpath, evpath+'/work/semblance_%s/' %kiter)
                    os.system('cp -r %s %s_%s/' % (evpath+'/*.config*', evpath+'/work/semblance', kiter))
                    os.system('cp -r %s %s_%s/' % (evpath+'/*.origin*', evpath+'/work/semblance', kiter))
                    os.system('cp -r %s %s_%s/' % (evpath+'/*.syn*', evpath+'/work/semblance', kiter))


def processLoop(traces=None, stations=None, cluster=None):
    force = False
    move = False
    for argv in sys.argv:
        if argv[-7:] == '--force':
            force = True
            evpath_emps = evpath_emp[:-7]
        elif argv[-12:] == '--force-move':
            force = False
            move = True
            evpath_emps = evpath_emp[:-12]
        else:
            evpath_emps = evpath_emp

    if force is False:
        check_is_empty(evpath, move=move)

    C = config.Config(evpath, eventpath_emp=evpath_emps)
    Origin = C.parseConfig('origin')
    C = config.Config(evpath, eventpath_emp=evpath_emps)

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
        Syn_in_emp = C.parseConfig('syn')
        syn_in_emp = SynthCfg(Syn_in)

    yaml_file = C.parseConfig('yaml')
    cfg = guts.load(filename=yaml_file[0])
    Config_cluster = C.parseConfig('config')
    cfgn = ConfigObj(dict=Config_cluster)
    phases = cfg.config.ttphases

    if cfg.config_data.pyrocko_download is True:
        Meta = C.readpyrockostations(phases)
        if len(Meta) == 0:
            Meta = C.readpyrockostations('P')
    elif cfg.config.colesseo_input is True:
        scenario = guts.load(filename=cfg.config.colosseo_scenario_yml())
        scenario_path = cfg.config.colosseo_scenario_yml[:-12]
        Meta = C.readcolosseostations(scenario_path)
    else:
        Meta = C.readMetaInfoFile()
    Folder = C.createFolder()
    C.writeConfig(Config_cluster, Origin, Folder)
    filter = cfg.config_filter
    if filter.forerun > 0:
        ntimes = int((filter.forerun + filter.duration) /
                     filter.step)
    else:
        ntimes = int((filter.duration) / filter.step)
    origin = OriginCfg(Origin)

    if cfg.config_data.colesseo_input is True:
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

    if cfg.config_weight.correct_shifts_empirical is True:

        Origin_emp = C.parseConfig('origin_emp')
        origin_emp = OriginCfg(Origin_emp)
        ev_emp = Event(origin_emp.lat(), origin_emp.lon(), origin_emp.depth(),
                       origin_emp.time(), strike=strike, dip=dip, rake=rake)
    filtername = filter.name[0]
    Logfile.add('filtername = ' + filtername)

    XDict = OrderedDict()
    RefDict = OrderedDict()
    SL = OrderedDict()
    refshifts_global = []
    newFreq = str(filter.newFrequency)
    xcorrnetworks = cfgn.String('networks').split(',')

    if cfg.config_xcorr.xcorr is True:

        fobjreferenceshiftname = newFreq + '_' + filtername + '.refpkl'
        rp = os.path.join(Folder['semb'], fobjreferenceshiftname)
        fobjreferenceshiftnameemp = newFreq + '_' + filtername + 'emp' + '.refpkl'
        rpe = os.path.join(Folder['semb'], fobjreferenceshiftnameemp) + '.shift'
        fobjpickleshiftname = newFreq + '_' + filtername + '.xcorrpkl'
        ps = os.path.join(Folder['semb'], fobjpickleshiftname)
        if cfg.config_syn.synthetic_test is False:
            if cfg.config_data.quantity == 'displacement':
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
                SL[i] = len(Config_cluster[i].split('|'))
        else:
            if cfg.config_syn.synthetic_test is False:
                if cfg.config_data.quantity == 'displacement':
                    try:
                        traces = io.load(evpath+'/data/traces_rotated.mseed')
                    except Exception:
                        traces = io.load(evpath+'/data/traces_restituted.mseed')
                else:
                    traces = io.load(evpath+'/data/traces_velocity.mseed')
            SL = {}
            for i in xcorrnetworks:
                W = {}
                network = cfgn.String(i).split('|')
                FilterMeta = ttt.filterStations(Meta, Config_cluster, Origin,
                                                network, cfg)
                arrayfolder = os.path.join(Folder['semb'], i)
                if os.access(arrayfolder, os.F_OK) is False:
                    os.makedirs(arrayfolder)
                if cfg.config_data.pyrocko_download is True:
                    # TODO check seperate xcoor nescessity
                    A = Xcorr(ev, FilterMeta, evpath, Config_cluster, Syn_in,
                              arrayfolder)
                print("run Xcorr")
                phase = phases[0]
                W, triggerobject = A.runXcorr(phase, traces, cfg)
                XDict[i] = W
                RefDict[i] = triggerobject.tdiff
                SL[i] = len(network)
                for j in range(0, len(FilterMeta)):
                    refshifts_global.append(triggerobject.tdiff)

    else:
        fobjreferenceshiftname = newFreq + '_' + filtername + '.refpkl'
        rp = os.path.join(Folder['semb'], fobjreferenceshiftname)
        fobjreferenceshiftnameemp = newFreq + '_' + filtername + 'emp' + '.refpkl'
        rpe = os.path.join(Folder['semb'], fobjreferenceshiftnameemp) + '.shift'
        fobjpickleshiftname = newFreq + '_' + filtername + '.xcorrpkl'
        ps = os.path.join(Folder['semb'], fobjpickleshiftname)

        refshift = 0
        SL = {}
        for i in xcorrnetworks:
            W = {}
            network = cfgn.config.String(i).split('|')
            FilterMeta = ttt.filterStations(Meta, Config_cluster,
                                            Origin, network)
            arrayfolder = os.path.join(Folder['semb'], i)
            if os.access(arrayfolder, os.F_OK) is False:
                os.makedirs(arrayfolder)
            XDict[i] = FilterMeta
            RefDict[i] = refshift
            SL[i] = len(network)
            for j in range(0, len(FilterMeta)):
                refshifts_global.append(refshift)
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
                        % (Config_cluster['networks']))
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
            Config_cluster['networks'] = nnl
            break

    for j in range(3, 0, -1):
        time.sleep(1)
        Logfile.red('Start processing in %d seconds ' % (j))

    wd = Origin['depth']
    start, stop, step = cfg.config_geometry.depths

    start = int(start)
    stop = int(stop)+1
    step_depth = int(step)
    filters = cfg.config_filter.filters
    filters = int(filters)
    Logfile.add('working on ' + Config_cluster['networks'])
    if cfg.config_weight.correct_shifts_empirical is True:
        emp_loop = True
    else:
        emp_loop = False

# ==================================loop over phases======================
    for phase in phases:
        if phase is 'P':
            desired = 'Z'
        if phase is 'S':
            desired = 'T'
            rpe = rpe +'_S'
        # ==================================loop over filter setups=====
        for filterindex in xrange(0, filters):
            # ==================================loop over depth=======
            for depthindex in xrange(start, stop, step_depth):

                workdepth = float(wd) + depthindex
                if cfg.config_geometry.dimz == 0:
                    Origin['depth'] = workdepth
                ev = Event(Origin['lat'], Origin['lon'], Origin['depth'],
                           Origin['time'], strike=strike, dip=dip, rake=rake)
                Logfile.add('WORKDEPTH: ' + str(Origin['depth']))
                networks = Config_cluster['networks'].split(',')

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

                    network = Config_cluster[i].split('|')

                    Logfile.add('network: ' + str(network))
                    FilterMeta = ttt.filterStations(Meta, Config, Origin,
                                                    network, cfg)

                    W = XDict[i]
                    refshift = RefDict[i]
                    for j in range(0, len(FilterMeta)):
                        if cfg.config_weight.correct_shifts is False:
                            refshift = refshift*0.
                        refshifts.append(refshift)

                    Logfile.add('BOUNDING BOX DIMX: %s  DIMY: %s  GRIDSPACING:\
                                %s \n' % (str(cfg.config_geometry.dimx),
                                          str(cfg.config_geometry.dimy),
                                          str(cfg.config_geometry.dimz)))

                    Logfile.red('Calculating Traveltime Grid')
                    t1 = time.time()

                    isParallel = False
                    TTTGridMap = []
                    mint = []
                    maxt = []
                    ttt_model = cfg.config.traveltime_model
                    try:
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

                        except:
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
                                if cfg.config_geometry.dimz != 0:
                                    t1 = time.time()
                                    ttt.calcTTTAdv_cube(Config, FilterMeta[i],
                                                        Origin, i, arrayname,
                                                        W, refshift, phase)

                                    Logfile.add('ttt.calcTTTAdv : '
                                                + str(time.time() - t1)
                                                + ' sec.')
                                else:
                                    t1 = time.time()
                                    ttt.calcTTTAdv(cfg, FilterMeta[i],
                                                   Origin,
                                                   i, arrayname, W, refshift,
                                                   phase)

                                    Logfile.add('ttt.calcTTTAdv : '
                                                + str(time.time() - t1)
                                                + ' sec.')
                        assert len(FilterMeta) > 0
                        if cfg.config_geometry.dimz != 0:
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

                    if cfg.config_weight.correct_shifts_empirical is True:
                        TTTGridMap_emp = []
                        Logfile.add('BOUNDING BOX DIMX: %s  DIMY: %s  GRIDSPACING:\
                                    %s \n' % (Config['dimx_emp'],
                                              Config['dimy_emp'],
                                              Config['gridspacing']))
                        try:
                            px = os.path.abspath(os.path.join(os.getcwd(),
                                                 os.pardir))
                            pdx = str('/tttgrid/tttgrid%s_%s_%s_%s_%s_emp.pkl'
                                      % (phase, ttt_model, ev_emp.time, arrayname,
                                         workdepth))
                            px_path = px+pdx
                            f = open(px_path, 'rb')

                            print("loading travel time tttgrid%s_%s_%s_%s_%s_emp.pkl"
                                  % (phase, ttt_model, ev.time, arrayname,
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
                                                   (Config, FilterMeta[i], Origin,
                                                    i, arrayname, W, refshift))

                                    po.close()
                                    po.join()
                            else:
                                for i in xrange(len(FilterMeta)):
                                    if cfg.config_geometry.dimz != 0:
                                        t1 = time.time()
                                        ttt.calcTTTAdv_cube(Config,
                                                            FilterMeta[i],
                                                            Origin_emp, i,
                                                            arrayname,
                                                            W, refshift, phase,
                                                            flag_rpe=True)

                                        Logfile.add('ttt.calcTTTAdv : '
                                                    + str(time.time() - t1)
                                                    + ' sec.')
                                    else:
                                        t1 = time.time()
                                        ttt.calcTTTAdv(Config, FilterMeta[i],
                                                       Origin_emp,
                                                       i, arrayname, W,
                                                       refshift,
                                                       phase, flag_rpe=True)

                                        Logfile.add('ttt.calcTTTAdv : '
                                                    + str(time.time() - t1)
                                                    + ' sec.')
                            assert len(FilterMeta) > 0
                            if cfg.config_geometry.dimz != 0:
                                TTTGridMap_emp = deserializer.deserializeTTT_cube(len(FilterMeta), flag_rpe=True)
                                mint_emp, maxt_emp = deserializer.deserializeMinTMaxT(len(FilterMeta), flag_rpe=True)
                            else:
                                TTTGridMap_emp = deserializer.deserializeTTT(len(FilterMeta), flag_rpe=True)
                                mint_emp, maxt_emp = deserializer.deserializeMinTMaxT(len(FilterMeta), flag_rpe=True)

                            px = os.path.abspath(os.path.join(os.getcwd(),
                                                 os.pardir))
                            pdx = str('/tttgrid/tttgrid%s_%s_%s_%s_%s_emp.pkl'
                                      % (phase, ttt_model, ev_emp.time, arrayname,
                                         workdepth))
                            px_path = px+pdx
                            f = open(px_path, 'wb')
                            pickle.dump([TTTGridMap_emp, mint_emp, maxt_emp], f)
                            f.close()
                        t2 = time.time()
                        Logfile.red('%s took %0.3f s' % ('TTT', (t2-t1)))

                    switch = filterindex
                    tw = times.calculateTimeWindows(mint, maxt, cfg,
                                                    ev, switch)
                    if cfg.config_weight.correct_shifts_empirical is True:
                        tw_emp = times.calculateTimeWindows(mint_emp, maxt_emp,
                                                            cfg, ev_emp,
                                                            switch)

                        if cfg.config_data.pyrocko_download is True:

                            if cfg.config_weight.correct_shifts_empirical_synthetic is True:
                                Wd_emp = waveform.readWaveformsPyrockodummy(FilterMeta,
                                                                            tw,
                                                                            evpath,
                                                                            ev)
                            elif cfg.config_data.quantity == 'displacement':
                                Wd_emp = waveform.readWaveformsPyrocko_restituted(
                                    FilterMeta, tw, evpath, ev_emp, desired)
                            else:
                                Wd_emp = waveform.readWaveformsPyrocko(FilterMeta,
                                                                       tw_emp,
                                                                       evpath_emps,
                                                                       ev_emp,
                                                                       desired)
                        elif cfg.config_data.colesseo_input is True:
                            Wd_emp = waveform.readWaveforms_colesseo(FilterMeta,
                                                                     tw_emp,
                                                                     evpath_emps,
                                                                     ev_emp, C)
                        else:
                            Wd_emp = waveform.readWaveforms(FilterMeta, tw_emp,
                                                            evpath_emps, ev_emp)
                        if cfg.config_weight.correct_shifts_empirical_synthetic is True\
                           or cfg.config_filter.dynamic_filter is True:
                            Wdf_emp = waveform.processdummyWaveforms(Wd_emp, Config,
                                                                 Folder, arrayname,
                                                                 FilterMeta, ev_emp,
                                                                 switch, W)
                            Wdfs_emp.extend(Wdf_emp)
                        else:
                            if switch == 0:
                                ff1 = cfg.config_filter.flo[switch]
                                ff2 = cfg.config_filter.fhi[switch]
                                switchs = "l0"
                            else:
                                f1 = str('filter.flo%s()'% str(filterindex+1))
                                ff1 = eval(f1)
                                f2 = str('filter.fhi%s()'% str(filterindex+1))
                                ff2 = eval(f2)
                                switchs = "h1"
                            ps_wdf_emp = os.path.join(Folder['semb'],
                                                      "fobjpickle_process_emp\
                                                      _%s_%s%s"
                                                      % (arrayname, ff1, ff2))
                            if cfg.config_data.load_wdf is True:
                                try:
                                    f = open(ps_wdf_emp, 'rb')
                                    Wdf_emp = pickle.load(f)
                                except Exception:
                                    Wdf_emp = waveform.processdummyWaveforms(Wd_emp,
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
                                Wdf_emp = waveform.processdummyWaveforms(Wd_emp,
                                                                    Config,
                                                                    Folder,
                                                                    arrayname,
                                                                    FilterMeta,
                                                                    ev_emp,
                                                                    switch, W)
                            Wdfs_emp.extend(Wdf_emp)
                    if cfg.config_data.pyrocko_download is True:

                        if cfg.config_data.quantity  == 'displacement':
                            Wd = waveform.readWaveformsPyrocko_restituted(
                                FilterMeta, tw, evpath, ev, desired)
                        elif cfg.config_syn.synthetic_test is True:
                            Wd = waveform.readWaveformsPyrockodummy(FilterMeta,
                                                                    tw, evpath,
                                                                    ev)
                        else:
                            Wd = waveform.readWaveformsPyrocko(FilterMeta, tw,
                                                               evpath, ev,
                                                               desired)
                    elif cfg.config_data.colesseo_input is True:
                        Wd = waveform.readWaveforms_colesseo(FilterMeta, tw,
                                                             evpath, ev, C)
                    else:
                        Wd = waveform.readWaveforms(FilterMeta, tw, evpath, ev)

                    if cfg.config_syn.synthetic_test is True\
                       or cfg.config_filter.dynamic_filter is True:
                        Wdf = waveform.processdummyWaveforms(Wd, Config,
                                                             Folder, arrayname,
                                                             FilterMeta, ev,
                                                             switch, W)
                        Wdfs.extend(Wdf)
                    else:
                        if switch == 0:
                            ff1 = cfg.config_filter.flo[switch]
                            ff2 = cfg.config_filter.fhi[switch]
                            switchs = "l0"
                        else:
                            f1 = str('filter.flo%s()'% str(filterindex+1))
                            ff1 = eval(f1)
                            f2 = str('filter.fhi%s()'% str(filterindex+1))
                            ff2 = eval(f2)
                            switchs = "h1"
                        ps_wdf = os.path.join(Folder['semb'], "fobjpickle_process_%s_%s%s" % (arrayname, ff1, ff2))
                        if cfg.config_data.load_wdf is True:
                            try:
                                f = open(ps_wdf, 'rb')
                                Wdf = pickle.load(f)
                            except:
                                Wdf = waveform.processdummyWaveforms(Wd, Config,
                                                                Folder,
                                                                arrayname,
                                                                FilterMeta,
                                                                ev, switch, W)

                                fobj_proc = open(ps_wdf, 'wb')
                                pickle.dump(Wdf, fobj_proc)
                                f = open(ps_wdf, 'rb')
                                Wdf = pickle.load(f)
                        else:
                            Wdf = waveform.processdummyWaveforms(Wd, Config, Folder,
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
                        step = cfg.config_filter.step
                    else:
                        s1 = str('cfg.config.step_f%s()')% str(filterindex+1)
                        step = eval(s1)

                    if cfg.config_filter.forerun > 0:
                        ntimes = int((cfg.config_filter.forerun +
                                      cfg.config_filter.duration) / step)
                    else:
                        ntimes = int((cfg.config_filter.duration) / step)
                    if cfg.config_weight.combine_all is False:

                        if cfg.config.optimize is True:
                            optim.solve(counter, Config, Wdf, FilterMeta,
                                        mint, maxt, TTTGridMap, Folder,
                                        Origin, ntimes, switch, ev,
                                        arrayfolder, syn_in, refshifts, phase,
                                        rpe+str(arrayname)+switchs, flag_rpe)
                        else:
                            if cfg.config_weight.correct_shifts_empirical is True:
                                if cfg.config_weight.correct_shifts_empirical_run is True:
                                    winlen_emp = cfg.config.winlen_emp()
                                    step_emp = cfg.config.step_emp()
                                    if cfg.config_filter.forerun > 0:
                                        ntimes_emp = int((cfg.config_filter.forerun_emp + cfg.config_filter.duration_emp)/step_emp)
                                    else:
                                        ntimes_emp = int((cfg.config_filter.duration_emp) / step_emp)
                                    f = open(os.path.abspath(os.path.join(os.getcwd(),
                                                             os.pardir))+'/tttgrid/tttgrid%s_%s_%s_%s_%s_emp.pkl'
                                             % (phase, ttt_model, ev_emp.time,
                                                arrayname, workdepth), 'rb')
                                    print("loading travel time grid%s_%s_%s_%s_%s_emp.pkl"
                                          % (phase, ttt_model, ev_emp.time, arrayname,
                                             workdepth))
                                    TTTGridMap_emp, mint_emp, maxt_emp = pickle.load(f)
                                    f.close()
                                    if cfg.config_data.pyrocko_download is True:

                                        if cfg.config_weight.correct_shifts_empirical_synthetic is True:
                                            Wd_emp = waveform.readWaveformsPyrockodummy(FilterMeta,
                                                                                    tw,
                                                                                    evpath,
                                                                                    ev)
                                        elif cfg.config_data.quantity  == 'displacement':
                                            Wd_emp = waveform.readWaveformsPyrocko_restituted(
                                                FilterMeta, tw, evpath, ev_emp, desired)
                                        else:
                                            Wd_emp = waveform.readWaveformsPyrocko(FilterMeta,
                                                                                   tw_emp,
                                                                                   evpath_emps,
                                                                                   ev_emp,
                                                                                   desired)
                                    elif cfg.config_data.colesseo_input is True:
                                        Wd_emp = waveform.readWaveforms_colesseo(FilterMeta,
                                                                                 tw_emp,
                                                                                 evpath_emps,
                                                                                 ev_emp, C)
                                    else:
                                        Wd_emp = waveform.readWaveforms(FilterMeta, tw_emp,
                                                                        evpath_emps, ev_emp)
                                    flag_rpe = True
                                    nstats = stations_per_array
                                    if switch == 0:
                                        switchs = "l0"
                                    else:
                                        switchs = "h1"

                                    arraySemb, weight, array_center = sembCalc.doCalc(
                                        counter, cfg, Wd_emp, FilterMeta,
                                        mint_emp, maxt_emp,
                                        TTTGridMap_emp, Folder, Origin_emp,
                                        ntimes_emp, switch, ev_emp,
                                        arrayfolder, syn_in_emp, refshifts,
                                        phase, rpe+str(arrayname)+switchs, flag_rpe,
                                        nstats)

                            if switch == 0:
                                switchs = "l0"
                            else:
                                switchs = "h1"

                            if cfg.config_weight.correct_shifts_empirical is True:

                                if sys.version_info.major >= 3:
                                    f = open(rpe+str(arrayname)+switchs, 'rb')
                                else:
                                    f = open(rpe+str(arrayname)+switchs)
                                RefDict_empirical = pickle.load(f)
                                refshifts = RefDict_empirical

                                for j in range(0, len(Wd_emp)):
                                    if cfg.config_weight.correct_shifts is False:
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
                            if switch == 0:
                                switchs = "l0"
                            else:
                                switchs = "h1"

                            if cfg.config_data.pyrocko_download is True:

                                if cfg.config_data.quantity  == 'displacement':
                                    Wd = waveform.readWaveformsPyrocko_restituted(
                                        FilterMeta, tw, evpath, ev, desired)
                                elif cfg.config_syn.synthetic_test is True:
                                    Wd = waveform.readWaveformsPyrockodummy(FilterMeta,
                                                                            tw, evpath,
                                                                            ev)
                                else:
                                    Wd = waveform.readWaveformsPyrocko(FilterMeta, tw,
                                                                       evpath, ev,
                                                                       desired)
                            elif cfg.config_data.colesseo_input is True:
                                Wd = waveform.readWaveforms_colesseo(FilterMeta, tw,
                                                                     evpath, ev, C)
                            else:
                                Wd = waveform.readWaveforms(FilterMeta, tw, evpath, ev)

                            try:
                                f = open(os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/tttgrid/tttgrid%s_%s_%s_%s_%s.pkl'
                                         % (phase, ttt_model, ev.time, arrayname,
                                            workdepth), 'rb')
                                TTTGridMap, mint, maxt = pickle.load(f)
                                f.close()
                                arraySemb, weight, array_center = sembCalc.doCalc(
                                    counter, cfg, Wd, FilterMeta, mintt, maxtt,
                                    TTTGridMap, Folder, Origin, ntimes, switch,
                                    ev, arrayfolder, syn_in, refshifts, phase,
                                    rpe+str(arrayname)+switchs, flag_rpe,
                                    len(FilterMeta))
                            except:
                                try:

                                    f = open(os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/tttgrid/tttgrid%s_%s_%s_%s_%s.pkl'
                                             % (phase, ttt_model, ev.time, arrayname,
                                                workdepth), 'rb')
                                    TTTGridMap, mint, maxt = pickle.load(f)
                                    f.close()
                                    arraySemb, weight, array_center = sembCalc.doCalc(
                                        counter, cfg, Wd, FilterMeta, mintt, maxtt,
                                        TTTGridMap, Folder, Origin, ntimes, switch,
                                        ev, arrayfolder, syn_in, refshifts, phase,
                                        rpe+str(arrayname)+switchs, flag_rpe,
                                        len(FilterMeta))
                                except:
                                    try:
                                        f = open(os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/tttgrid/tttgrid%s_%s_%s_%s_%s.pkl'
                                                 % (phase, ttt_model, ev.time, arrayname,
                                                    workdepth), 'rb')
                                        TTTGridMap, mint, maxt = pickle.load(f)
                                        f.close()
                                        arraySemb, weight, array_center = sembCalc.doCalc(
                                            counter, cfg, Wd, FilterMeta, mintt, maxtt,
                                            TTTGridMap, Folder, Origin, ntimes, switch,
                                            ev, arrayfolder, syn_in, refshifts, phase,
                                            rpe+str(arrayname)+switchs, flag_rpe,
                                            len(FilterMeta))
                                    except:
                                        try:
                                            TTTGridMap = []
                                            f = open(os.path.abspath(os.path.join(os.getcwd(),
                                                     os.pardir))+'/tttgrid/tttgrid%s_%s_%s_%s_%s.pkl'
                                                     % (phase, ttt_model, ev.time, arrayname,
                                                        workdepth), 'rb')
                                            print("loading travel time grid%s_%s_%s_%s_%s.pkl"
                                                  % (phase, ttt_model, ev.time, arrayname,
                                                      workdepth))
                                            TTTGridMap, mintt, maxtt = pickle.load(f)
                                            f.close()
                                            arraySemb, weight, array_center = sembCalc.doCalc(
                                                counter, cfg, Wd, FilterMeta, mintt, maxtt,
                                                TTTGridMap, Folder, Origin, ntimes, switch,
                                                ev, arrayfolder, syn_in, refshifts, phase,
                                                rpe+str(arrayname)+switchs, flag_rpe,
                                                len(FilterMeta))
                                        except:
                                            f = open(os.path.abspath(os.path.join(os.getcwd(),
                                                     os.pardir))+'/tttgrid/tttgrid%s_%s_%s_%s_%s.pkl'
                                                     % (phase, ttt_model, ev.time, arrayname,
                                                        workdepth), 'rb')
                                            print("loading travel time grid%s_%s_%s_%s_%s.pkl"
                                                  % (phase, ttt_model, ev.time, arrayname,
                                                      workdepth))
                                            TTTGridMap, mintt, maxtt = pickle.load(f)
                                            f.close()
                                        arraySemb, weight, array_center = sembCalc.doCalc(
                                            counter, cfg, Wd, FilterMeta, mintt, maxtt,
                                            TTTGridMap, Folder, Origin, ntimes, switch,
                                            ev, arrayfolder, syn_in, refshifts, phase,
                                            rpe+str(arrayname)+switchs, flag_rpe,
                                            len(FilterMeta))
                            weights.append(weight)
                            array_centers.append(array_center)
                            ASL.append(arraySemb)
                            sembCalc.writeSembMatricesSingleArray(arraySemb,
                                                                  cfg,
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
                    stations_per_array.append(len(Wd))
                    TTTgrids.update(TTTGridMap)
                    mints.append(mint)
                    maxts.append(maxt)
                    FilterMetas[len(FilterMetas):] = FilterMeta
                    TTTGridMap = []

                if cfg.config_weight.combine_all is True:
                    if cfg.config_data.pyrocko_download is True:
                        if cfg.config_syn.synthetic_test is True:
                            Wd = waveform.readWaveformsPyrockodummy(
                                    FilterMetas, tw, evpath, ev)
                        else:
                            if cfg.config_data.quantity  == 'displacement':
                                Wd = waveform.readWaveformsPyrocko_restituted(
                                    FilterMetas, tw, evpath, ev, desired)
                            else:
                                Wd = waveform.readWaveformsPyrocko(FilterMetas,
                                                                   tw, evpath,
                                                                   ev, desired)
                                print(Wd, "Wd")
                    elif cfg.config_data.colesseo_input is True:
                        Wd = waveform.readWaveforms_colesseo(FilterMetas, tw,
                                                             evpath, ev, C)
                    else:
                        Wd = waveform.readWaveforms(FilterMetas, tw, evpath,
                                                    ev)
                    if cfg.config_syn.synthetic_test is True:
                        Wdf = waveform.processdummyWaveforms(Wd, Config,
                                                             Folder, arrayname,
                                                             FilterMetas, ev,
                                                             switch, W)
                    else:
                        if switch == 0:
                            ff1 = cfg.config_filter.flo[switch]
                            ff2 = cfg.config_filter.flo[switch]
                            switchs = "l0"
                        else:
                            f1 = str('filter.flo%s()'% str(filterindex+1))
                            ff1 = eval(f1)
                            f2 = str('filter.fhi%s()'% str(filterindex+1))
                            ff2 = eval(f2)
                            switchs = "h1"
                        ps_wdf = os.path.join(Folder['semb'],
                                              "fobjpickle_process_%s_%s%s_\
                                              combined" % (arrayname, ff1, ff2))
                        if cfg.config_data.load_wdf is True:
                            try:
                                f = open(ps_wdf, 'rb')
                                Wdf = pickle.load(f)
                                print('loaded wdf')
                            except:
                                Wdf = waveform.processdummyWaveforms(Wd, Config,
                                                                Folder,
                                                                arrayname,
                                                                FilterMetas,
                                                                ev, switch, W)

                                fobj_proc = open(ps_wdf, 'wb')
                                pickle.dump(Wdf, fobj_proc)
                                print('dumped wdf')
                        else:
                            Wdf = waveform.processdummyWaveforms(Wd, Config, Folder,
                                                            arrayname,
                                                            FilterMetas,
                                                            ev, switch, W)

                    mint = num.min(mints)
                    maxt = num.max(maxts)
                    nstats = stations_per_array
                    flag_rpe = False
                    if switch == 0:
                        ff1 = cfg.config_filter.flo[switch]
                        ff2 = cfg.config_filter.flo[switch]
                        switchs = "l0"
                    else:
                        f1 = str('filter.flo%s()'% str(filterindex+1))
                        ff1 = eval(f1)
                        f2 = str('ff2 = filter.fhi%s()'% str(filterindex+1))
                        ff2 = eval(f2)
                        switchs = "h1"

                    if cfg.config_weight.bootstrap_array_weights is False:
                        print(Wd, "Wd2")

                        arraySemb, weight, array_center = sembCalc.doCalc(
                            counter, cfg, Wd, FilterMetas, mint, maxt,
                            TTTgrids, Folder, Origin, ntimes, switch,
                            ev, arrayfolder, syn_in, refshifts_global, phase,
                            rpe+str(arrayname)+switchs, flag_rpe, nstats)
                        ASL.append(arraySemb)
                        weights.append(weight)
                        array_centers.append(array_center)
                        sembCalc.writeSembMatricesSingleArray(arraySemb,
                                                              Config, Origin,
                                                              arrayfolder,
                                                              ntimes, switch,
                                                              phase)
                    else:
                        nboot = cfg.config_weight.n_bootstrap
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
                            array_index = 0
                            boot_shifts = num.zeros(num.shape(SL))
                            for i in xcorrnetworks:
                                len_network = SL[i]
                                array_shift = rstate.uniform(-1*cfg.config_weight.shift_max,
                                                 cfg.config_weight.shift_max)
                                boot_shifts[array_index:array_index+len_network] = array_shift
                                array_index = array_index+len_network

                            arraySemb, weight, array_center = sembCalc.doCalc(
                                counter, cfg, Wdf, FilterMetas, mint, maxt,
                                TTTgrids, Folder, Origin, ntimes, switch,
                                ev, arrayfolder, syn_in, refshifts_global,
                                phase, rpe+str(arrayname)+switchs, flag_rpe, nstats,
                                bs_weights=ws, boot_shifts=boot_shifts)

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
                                # Todo: Add optim for combined semblance
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
                if cfg.config.optimize_all is True:
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
                    sembmax, tmp = sembCalc.collectSemb(ASL, cfg, origin,
                                                        Folder,
                                                        ntimes, len(networks),
                                                        switch, array_centers,
                                                        phase,
                                                        time=Origin['time'],
                                                        Origin=Origin)
                    if cfg.config_weight.weight_by_noise is True:
                        sembCalc.collectSembweighted(ASL, cfg, Origin,
                                                     Folder, ntimes,
                                                     len(networks), switch,
                                                     weights,
                                                     time=Origin['time'])

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
