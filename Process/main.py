import os
import sys
from pyrocko import guts
import logging
import shutil
import time
import multiprocessing
from optparse import OptionParser
if sys.version_info.major >= 3:
    import _pickle as pickle
    xrange = range
else:
    import cPickle as pickle
from pyrocko import obspy_compat

sys.path.append('../Common/')
import optim
import Basic
import Globals
import Logfile
from Program import MainObj
from ConfigFile import ConfigObj, FilterCfg, OriginCfg, SynthCfg
from collections import OrderedDict

sys.path.append('../tools/')
import config
from config import Event
import deserializer
import ttt
import sembCalc
import waveform
import times
from array_crosscorrelation_v4  import Xcorr, cmpFilterMetavsXCORR
import numpy as num

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

    global  evpath

    parser = OptionParser(usage="%prog -f Eventpath ")
    parser.add_option("-f", "--evpath", type="string", dest="evpath", help="evpath")

    (options, args) = parser.parse_args()

    if options.evpath == None:
        parser.error("non existing eventpath")
        return False

    evpath = options.evpath
    Globals.setEventDir(evpath)
    return True


def processLoop():

    #==================================get meta info========================
    C = config.Config(evpath)
    Origin = C.parseConfig('origin')
    try:
        Syn_in = C.parseConfig('syn')
        syn_in = SynthCfg(Syn_in)
    except:
        pass
    Config = C.parseConfig('config')

    cfg = ConfigObj(dict=Config)
    phases = cfg.Str('ttphases')
    phases = phases.split(',')

    if cfg.pyrocko_download() == True:
        Meta = C.readpyrockostations()
    elif cfg.colesseo_input() == True:
        scenario = guts.load(filename=cfg.colosseo_scenario_yml())
        scenario_path = cfg.colosseo_scenario_yml()[:-12]
        Meta = C.readcolosseostations(scenario_path)
    else:
        Meta = C.readMetaInfoFile()
    Folder = C.createFolder()
    C.writeConfig(Config,Origin,Folder)

    filter = FilterCfg(Config)
    if cfg.UInt('forerun')>0:
        ntimes = int((cfg.UInt('forerun') + cfg.UInt('duration') ) / cfg.UInt('step') )
    else:
        ntimes = int((cfg.UInt('duration') ) / cfg.UInt('step') )
    origin = OriginCfg(Origin)

    if cfg.colesseo_input() == True:
        from pyrocko import util
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

    filtername = filter.filterName()
    Logfile.add('filtername = ' + filtername)

    XDict = {}
    RefDict = {}
    SL = {}
    if cfg.Int('xcorr') == 1:

        newFreq = str(filter.newFrequency())
        fobjreferenceshiftname = newFreq + '_' + filtername + '.refpkl'
        rp= os.path.join(Folder['semb'], fobjreferenceshiftname)
        fobjpickleshiftname= newFreq + '_' + filtername + '.xcorrpkl'
        ps= os.path.join(Folder['semb'], fobjpickleshiftname)

        if(os.path.isfile(rp) and os.path.getsize(rp) != 0 and os.path.isfile(ps) and os.path.getsize(ps) != 0):
            Logfile.add('file exits : ' + rp)
            Logfile.add('load refshifts')

            if sys.version_info.major >= 3:
                f= open(rp, 'rb')
            else:
                f= open(rp)

            RefDict = pickle.load(f)
            if sys.version_info.major >= 3:
                x= open(ps, 'rb')
            else:
                x= open(ps)
            XDict= pickle.load(x)
            xcorrnetworks = cfg.String('networks').split(',')

            for i in xcorrnetworks:
                SL[i] = len(Config[i].split('|'))
        else:
            SL = {}
            xcorrnetworks = cfg.String('networks').split(',')
            for i in xcorrnetworks:
                W = {}
                refshift= 0
                network = cfg.String(i).split('|')
                FilterMeta = ttt.filterStations(Meta,Config,Origin,network)
                arrayfolder = os.path.join(Folder['semb'],i)

                if os.access(arrayfolder,os.F_OK) == False:
                   os.makedirs(arrayfolder)
                if cfg.pyrocko_download() == True:
                    A = Xcorr(ev,FilterMeta,evpath,Config,Syn_in,arrayfolder)
                else:
                    A = Xcorr(ev,FilterMeta,evpath,Config,Syn_in,arrayfolder)

                print("run Xcorr")
                phase = phases[0]
                W,triggerobject= A.runXcorr(phase)

                XDict[i] = W
                RefDict[i] = triggerobject.tdiff
                SL[i] = len(network)

            if sys.version_info.major >= 3:
                fobjrefshift = open(rp,'wb')
            else:
                fobjrefshift = open(rp,'w')
            pickle.dump(RefDict,fobjrefshift)
            fobjrefshift.close()

            if sys.version_info.major >= 3:
                output = open(ps, 'wb')
            else:
                output = open(ps, 'w')
            pickle.dump(XDict, output)
            output.close()
        if sys.version_info.major >= 3:
            for i in sorted(XDict.keys()) :
                Logfile.red('Array %s has %3d of %3d Stations left' %(i,len(XDict[i]),SL[i]))
        else:
            for i in sorted(XDict.iterkeys()) :
                Logfile.red('Array %s has %3d of %3d Stations left' %(i,len(XDict[i]),SL[i]))

        logger.info('\033[31mFor proceeding without changes press enter or give new comma seperatet network list or quit for exit\033[0m')

        while True:
            if sys.version_info.major >= 3:
                nnl = input("please enter your choice: ")
            else:
                nnl = raw_input("please enter your choice: ")

            if len(nnl) == 0:
                if not Basic.question('Process all networks ?') : continue

                Logfile.red('This networks will be used for processing: %s' %(Config['networks']))
                break

            elif str(nnl) == 'quit':
                sys.exit()

            elif str(nnl) == 'rerun':
                event = os.path.join(*evpath.split('/')[-1:])

                try:
                   os.remove(rp)
                   os.remove(ps)

                except : pass

                mainfolder = os.path.join(os.path.sep,*evpath.split('/')[:-2])
                os.chdir(mainfolder)

                cmd =('%s arraytool.py process %s') %(sys.executable,event)
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
                    continue   # Illegal network : input again

                Logfile.add('This networks will be used for processing: %s' %(nnl))
                Config['networks'] = nnl
                break

        for i in range(3, 0, -1):
            time.sleep(1)
            Logfile.red('Start processing in %d seconds ' %(i))


    wd = Origin['depth']
    start,stop,step = cfg.String('depths').split(',')

    start = int(start)
    stop = int(stop)+1
    step = int(step)
    filters = cfg.String('filters')
    filters = int(filters)
    Logfile.add('working on ' + Config['networks'])

# ==================================loop over phases======================
    for phase in phases:
        if phase is 'P':
            desired = 'Z'
        if phase is 'S':
            desired = 'T'
        # ==================================loop over filter setups=====
        for filterindex in xrange(0, filters):
            # ==================================loop over depth=======
            for depthindex in xrange(start, stop, step):

                workdepth = float(wd) + depthindex
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
                FilterMetas = []
                TTTgrids = OrderedDict()
                mints = []
                maxts = []

                for i in networks:

                    arrayname = i
                    arrayfolder = os.path.join(Folder['semb'], arrayname)

                    network = Config[i].split('|')
                    Logfile.add('network: ' + str(network))

                    FilterMeta = ttt.filterStations(Meta, Config, Origin,
                                                    network)

                    if len(FilterMeta) < 3:
                        continue

                    W = XDict[i]
                    refshift = RefDict[i]

                    FilterMeta = cmpFilterMetavsXCORR(W, FilterMeta)

                    Logfile.add('BOUNDING BOX DIMX: %s  DIMY: %s  GRIDSPACING: %s \n'
                             %(Config['dimx'],Config['dimy'],Config['gridspacing']))

                    Logfile.red('Calculating Traveltime Grid')
                    t1 = time.time()

                    isParallel = False
                    TTTGridMap = []
                    mint = []
                    maxt = []
                    ttt_model = cfg.Str('traveltime_model')

                    try:
                        f = open('../tttgrid/tttgrid%s_%s_%s_%s_%s.pkl' %(phase, ttt_model,ev.time, arrayname, workdepth), 'rb')
                        print("loading travel time grid%s_%s_%s_%s_%s.pkl" %(phase, ttt_model, ev.time, arrayname, workdepth))
                        TTTGridMap,mint,maxt = pickle.load(f)
                        f.close()
                        print("loading of travel time grid sucessful")
                    except:
                        print("loading of travel time grid unsucessful, will now calculate the grid:")
                        if isParallel :
                           maxp = 6
                           po = multiprocessing.Pool(maxp)

                           for i in xrange(len(FilterMeta)):
                               po.apply_async(ttt.calcTTTAdv,(Config,FilterMeta[i],Origin,i,arrayname,W,refshift))

                               po.close()
                               po.join()
                        else:
                            for i in xrange(len(FilterMeta)):
                                t1 = time.time()
                                ttt.calcTTTAdv(Config,FilterMeta[i],Origin,i,arrayname,W,refshift,phase)

                                Logfile.add('ttt.calcTTTAdv : ' + str(time.time() - t1) + ' sec.')

                        assert len(FilterMeta) > 0

                        TTTGridMap = deserializer.deserializeTTT(len(FilterMeta))
                        mint,maxt = deserializer.deserializeMinTMaxT(len(FilterMeta))
                        f = open('../tttgrid/tttgrid%s_%s_%s_%s_%s.pkl' %(phase, ttt_model,ev.time, arrayname, workdepth), 'wb')
                        print("dumping the traveltime grid for this array")
                        pickle.dump([TTTGridMap, mint, maxt], f)
                        f.close()

                    t2 = time.time()
                    Logfile.red('%s took %0.3f s' % ('TTT', (t2-t1)))

                    switch = filterindex

                    tw = times.calculateTimeWindows(mint, maxt, Config,
                                                    ev, switch)
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
                        Wdf = waveform.processWaveforms(Wd, Config, Folder,
                                                        arrayname, FilterMeta,
                                                        ev, switch, W)
                        Wdfs.extend(Wdf)

                    C.writeStationFile(FilterMeta, Folder, counter)
                    Logfile.red('%d Streams added for Processing' % (len(Wd)))

                    t1 = time.time()
                    f = open('../tttgrid/tttgrid%s_%s_%s_%s_%s.pkl' % (phase,
                             ttt_model, ev.time, arrayname, workdepth), 'rb')
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
                            optim.solve(counter, Config, Wdf, FilterMeta, mint,
                                        maxt, TTTGridMap, Folder, Origin,
                                        ntimes, switch, ev, arrayfolder,
                                        syn_in)
                        else:

                            arraySemb, weight, array_center = sembCalc.doCalc(
                                counter, Config, Wdf, FilterMeta, mint, maxt,
                                TTTGridMap, Folder, Origin, ntimes, switch, ev,
                                arrayfolder, syn_in)
                            weights.append(weight)
                            array_centers.append(array_center)
                            ASL.append(arraySemb)
                            sembCalc.writeSembMatricesSingleArray(arraySemb,
                                                                  Config,
                                                                  Origin,
                                                                  arrayfolder,
                                                                  ntimes,
                                                                  switch)

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
                                    FilterMetas, tw, evpath, ev)
                            else:
                                Wd = waveform.readWaveformsPyrocko(FilterMetas,
                                                                   tw, evpath,
                                                                   ev)
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
                        Wdf = waveform.processWaveforms(Wd, Config, Folder,
                                                        arrayname, FilterMetas,
                                                        ev, switch, W)
                    mint = num.min(mints)
                    maxt = num.max(maxts)
                    if cfg.Bool('bootstrap_array_weights') is False:

                        arraySemb, weight, array_center = sembCalc.doCalc(
                            counter, Config, Wdf, FilterMetas, mint, maxt,
                            TTTgrids, Folder, Origin, ntimes, switch,
                            ev, arrayfolder, syn_in)
                        ASL.append(arraySemb)
                        weights.append(weight)
                        array_centers.append(array_center)
                        sembCalc.writeSembMatricesSingleArray(arraySemb,
                                                              Config, Origin,
                                                              arrayfolder,
                                                              ntimes, switch)
                    else:
                        nboot = cfg.Int('n_bootstrap')
                        ws = num.zeros((nboot, counter-1))
                        for ibootstrap in range(nboot):
                            f = rstate.uniform(0., 1., size=counter-1)
                            f[0] = 0.
                            f[-1] = 1.
                            f = num.sort(f)
                            g = f[1:] - f[:-1]
                            ws[ibootstrap, :] = g
                            arraySemb, weight, array_center = sembCalc.doCalc(
                                counter, Config, Wdf, FilterMetas, mint, maxt,
                                TTTgrids, Folder, Origin, ntimes, switch,
                                ev, arrayfolder, syn_in, bs_weights=g)

                            ASL.append(arraySemb)
                            weights.append(weight)
                            array_centers.append(array_center)
                            sembCalc.writeSembMatricesSingleArray(arraySemb,
                                                                  Config,
                                                                  Origin,
                                                                  arrayfolder,
                                                                  ntimes,
                                                                  switch)

                if cfg.optimize_all() is True:
                    import optim_csemb
                    sembmax = sembCalc.collectSemb(ASL,Config,Origin,Folder,ntimes,len(networks),switch)
                    optim_csemb.solve(counter, Config,Wdf,FilterMeta,mint,maxt,TTTGridMap,
                                                 Folder,Origin,ntimes,switch, ev,arrayfolder,
                                                 syn_in, ASL, sembmax, evpath, XDict, RefDict,
                                                 workdepth, filterindex, Wdfs)

                if ASL:
                    Logfile.red('collect semblance matrices from all arrays')
                    sembmax = sembCalc.collectSemb(ASL, Config, Origin, Folder,
                                                   ntimes, len(networks),
                                                   switch, array_centers)
                    if cfg.Bool('weight_by_noise') is True:
                        sembCalc.collectSembweighted(ASL, Config, Origin,
                                                     Folder, ntimes,
                                                     len(networks), switch,
                                                     weights)

    else:
        Logfile.red('Nothing to do  -> Finish')
    print("last work depth:")
    print(workdepth)
# ----------------------------------------------------------------------------


class ProcessMain(MainObj):

    def __init__(self) :
        initModule()

        MainObj.__init__(self, self, '0.3', 'process_run.log', 'process.log')

    # ---------------------------------------------------------------------------------------------

    def init(self) :
        file = 'ak135.model'

        if not os.path.isfile(file) :
           source = os.path.join('..','tools', file)
           Basic.checkFileExists(source, isAbort = True)
           shutil.copy(source, file)

        return True

    # --------------------------------------------------------------------------------------------

    def process(self) :
        import cProfile
        cProfile.run('processLoop()', filename='test.profile')
        return True

    def finish(self) :    pass

# -------------------------------------------------------------------------------------------------

def MainProc() :

    mainObj = ProcessMain()

    mainObj.run()

# -------------------------------------------------------------------------------------------------

isClient = False

if __name__ == "__main__":

   if not isClient :
      MainProc()
