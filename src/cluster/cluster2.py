import logging
from optparse import OptionParser
from pyrocko import guts
import os
import sys
import fnmatch
import random
import time
import shutil
import math

from palantiri.common import Basic
from palantiri.common import Globals
from palantiri.common import Logfile
from palantiri.common import DataTypes
from palantiri.common.Program import MainObj
from palantiri.common.ObspyFkt import loc2degrees
from palantiri.common.ConfigFile import ConfigObj
from pyrocko import model
from palantiri.tools import config
sys.setrecursionlimit(1500)

logger = logging.getLogger(sys.argv[0])
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter("%(message)s")

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)


class Station(object):

    def __init__(self, net, sta, loc, comp, lat=0, lon=0, ele=0, dip=0, azi=0,
                 gain=0, member=-1):

        self.net = net
        self.sta = sta
        self.loc = loc
        self.comp = comp
        self.lat = lat
        self.lon = lon
        self.ele = ele
        self.dip = dip
        self.azi = azi
        self.gain = gain
        self.member = member

    def getName(self):
        return self.net+'.'+self.sta+'.'+self.loc+'.'+self.comp

    def __str__(self):
        return('%s.%s.%s.%s') % (self.net, self.sta, self.loc, self.comp)

    def __eq__(self, other):
        return self.getName() == other.getName()


class Centroid(object):

    def __init__(self, lat, lon, rank):

        self.lat = lat
        self.lon = lon
        self.rank = rank


def readMetaInfoFile(EventPath):

    Logfile.red('Parsing MetaInfoFile')

    try:
        for i in os.listdir(EventPath):
            if fnmatch.fnmatch(i, 'metainfo-*.meta'):
                evfile = os.path.join(EventPath, i)

        MetaL = []

        Logfile.add(evfile)
        fobj = open(evfile, 'r')

        for i in fobj:
            line = i.split()

            net = line[0]
            sta = line[1]
            loc = line[2]
            comp = line[3]
            lat = line[4]
            lon = line[5]
            ele = line[6]
            dip = line[7]
            azi = line[8]
            gain = line[9]

            if fnmatch.fnmatch(comp, 'Z'):
                MetaL.append(Station(net, sta, loc, comp, lat, lon, ele, dip,
                                     azi, gain))

        Logfile.red('%d ENTRIES IN METAFILE FOUND' % (len(MetaL)))
    except Exception:
        Logfile.red('METAFILE NOT READABLE')

    return MetaL


def readpyrockostations(path, disp, cfg):

    if disp is True:
        stations = model.load_stations(path+'/data/stations_cluster.txt')
    else:
        stations = model.load_stations(path+'/data/stations_cluster.txt')
    MetaL = []
    phases = cfg.Str('ttphases')
    phases = phases.split(',')
    for phase in phases:
        if phase is 'P':
            desired = 'Z'
        if phase is 'S':
            desired = 'T'
    for sl in stations:
            for channel in sl.channels:
                if channel.name == desired:
                    MetaL.append(Station(str(sl.network), str(sl.station),
                                         str(sl.location), str(channel.name),
                                         str(sl.lat), str(sl.lon),
                                         str(sl.elevation), str(channel.dip),
                                         str(channel.azimuth),
                                         str(channel.gain)))
    return MetaL


def readcolosseostations(scenario_path):
    stations = model.load_stations(scenario_path+'/meta/stations.txt')
    MetaL = []
    for sl in stations:
            channel = sl.channels[2]
            MetaL.append(Station(str(sl.network), str(sl.station),
                                 str(sl.location), str(sl.channels[2])[:3],
                                 str(sl.lat), str(sl.lon),
                                 str(sl.elevation), str(channel.dip),
                                 str(channel.azimuth),
                                 str(channel.gain)))
    return MetaL


def createFolder(EventPath):

    Folder = {}
    Logfile.red('Create working environment')

    Folder['cluster'] = os.path.join(EventPath, 'cluster')

    if os.access(Folder['cluster'], os.F_OK) is False:
        os.makedirs(Folder['cluster'])

    if os.access(os.getcwd(), os.W_OK):
        basedir = os.path.join(EventPath, 'work')
        sembdir = os.path.join(basedir, 'semblance')
        ascdir = os.path.join(basedir, 'asc')
        mseeddir = os.path.join(basedir, 'mseed')

        Folder['base'] = basedir
        Folder['semb'] = sembdir
        Folder['asc'] = ascdir
        Folder['mseed'] = mseeddir
        Folder['event'] = EventPath

    else:
        Logfile.abort('create Folder: No write permissions for ' + os.getcwd())

    Folder['config'] = os.path.join(os.getcwd(), 'skeleton')
    return Folder


def filterStations(StationList, Config, Origin):
    F = []
    cfg = ConfigObj(dict=Config)

    minDist, maxDist = cfg.FloatRange('mindist', 'maxdist')
    origin = DataTypes.dictToLocation(Origin)

    Logfile.red('Filter stations with configured parameters')

    for i in StationList:
        sdelta = loc2degrees(origin, i)

        if sdelta > minDist and sdelta < maxDist:
            F.append(Station(i.net, i.sta, i.loc, i.comp, i.lat, i.lon, i.ele,
                             i.dip, i.azi, i.gain))
    Logfile.red('%d STATIONS LEFT IN LIST' % len(F))
    return F


def checkStationAroundInitialCentroid(station, Config, StationMetaList):

    cfg = ConfigObj(dict=Config)
    initDist = cfg.Float('initialstationdistance')
    counter = 0

    for i in StationMetaList:
        sdelta = loc2degrees(station, i)
        if sdelta < initDist:
            counter +=1

    return counter


def addOK(station, stationList, Config, MetaList):

    cfg = ConfigObj(dict=Config)
    minDist = 0
    minAround = cfg.UInt('minstationaroundinitialcluster')
    t = 0

    for i in stationList:
        sdelta = loc2degrees(station, i)

        if sdelta > minDist:
            aroundcounter = checkStationAroundInitialCentroid(station,
                                                              Config,
                                                              MetaList)

            if aroundcounter >= minAround:
                t = 1
            else:
                t = 0
                return t
        else:
            t = 0
            return t

    return t


def DeterminedInitialCentroids(Config, StationMetaList):
    Logfile.red('Begin initial centroid search')
    cfg = ConfigObj(dict=Config)

    initialCentroids = []
    usedIndexes = []
    random.seed(time.clock())

    if len(StationMetaList) == 0:
        Logfile.red('Empty station list')
        return initialCentroids

    MAX_TIME_ALLOWED = 350
    start = time.time()
    radi = 0
    if int(Config['maxcluster']) == 0:
        to = len(StationMetaList)-1
    else:
        to = int(Config['maxcluster'])
    while len(initialCentroids) < to:
            dist_centroids = float(Config['centroidmindistance'])

            randomIndex = radi = radi+random.randint(0, 5)
            if randomIndex > len(StationMetaList):
                randomIndex = random.randint(0, 5)
            redraw = True
            while redraw is True:
                if randomIndex in usedIndexes:
                    randomIndex = random.randint(0, len(StationMetaList)-1)
                else:
                    if len(usedIndexes) > 2:
                        for rdx in usedIndexes:
                            s1 = StationMetaList[randomIndex]
                            s2 = StationMetaList[rdx]
                            delta = loc2degrees(s1, s2)
                            if delta >= dist_centroids:
                                redraw = False
                    else:
                        redraw = False
            usedIndexes.append(randomIndex)

            around = checkStationAroundInitialCentroid(StationMetaList[randomIndex],
                                                       Config, StationMetaList)
            found = False

            if len(initialCentroids) == 0:
                initialCentroids.append(StationMetaList[randomIndex])
                found = True
                start = time.time()

            else:
                t = addOK(StationMetaList[randomIndex], initialCentroids,
                          Config, StationMetaList)
                if (time.time() - start) > MAX_TIME_ALLOWED:
                    break
                if t == 1:
                    if len(usedIndexes) > 1:
                        for rdx in usedIndexes:
                            s1=StationMetaList[randomIndex]
                            s2=StationMetaList[rdx]
                            delta=loc2degrees(s1, s2)
                            if delta >= dist_centroids:
                                initialCentroids.append(StationMetaList[randomIndex])
                                found = True
                    else:
                        initialCentroids.append(StationMetaList[randomIndex])
                        found = True

                else:
                    continue

            if found:
                initDist = cfg.Float('initialstationdistance')
                Logfile.red('found initial cluster %d' %(len(initialCentroids)))
                Logfile.red('centroid %s with %d stations around %s deegree' % (StationMetaList[randomIndex], around, initDist))

    print(len(initialCentroids))
    Logfile.red('Initial centroid search finished')
    return initialCentroids


def createRandomInitialCentroids(Config, StationMetaList):
    Logfile.red('Begin initial centroid search')
    cfg = ConfigObj(dict=Config)

    initialCentroids = []
    usedIndexes = []
    random.seed(time.clock())

    if len(StationMetaList) == 0:
       Logfile.red('Empty station list')
       return initialCentroids

    MAX_TIME_ALLOWED = 350
    start = time.time()
    if int(Config['maxcluster']) == 0:
        to = len(StationMetaList)-1
    else:
        to = int(Config['maxcluster'])
    while len(initialCentroids) < to:
            dist_centroids = float(Config['centroidmindistance'])

            randomIndex = random.randint(0, len(StationMetaList)-1)
            redraw = True
            while redraw is True:
                if randomIndex in usedIndexes:
                    randomIndex = random.randint(0, len(StationMetaList)-1)
                else:
                    if len(usedIndexes) > 2:
                        for rdx in usedIndexes:
                            s1 = StationMetaList[randomIndex]
                            s2 = StationMetaList[rdx]
                            delta = loc2degrees(s1, s2)
                            if delta >= dist_centroids:
                                redraw = False
                    else:
                        redraw = False
            usedIndexes.append(randomIndex)

            around = checkStationAroundInitialCentroid(StationMetaList[randomIndex],
                                                       Config, StationMetaList)
            found = False

            if len(initialCentroids) == 0:
                initialCentroids.append(StationMetaList[randomIndex])
                found = True
                start = time.time()

            else:
                t = addOK(StationMetaList[randomIndex], initialCentroids,
                          Config, StationMetaList)
                if (time.time() - start) > MAX_TIME_ALLOWED:
                    break
                if t == 1:
                    if len(usedIndexes) > 1:
                        for rdx in usedIndexes:
                            s1=StationMetaList[randomIndex]
                            s2=StationMetaList[rdx]
                            delta=loc2degrees(s1, s2)
                            if delta >= dist_centroids:
                                initialCentroids.append(StationMetaList[randomIndex])
                                found = True
                    else:
                        initialCentroids.append(StationMetaList[randomIndex])
                        found = True

                else:
                    continue

            if found:
                initDist = cfg.Float('initialstationdistance')
                Logfile.red('found initial cluster %d' % (len(initialCentroids)))
                Logfile.red('centroid %s with %d stations around %s deegree' % (StationMetaList[randomIndex], around, initDist))


    Logfile.red('Initial centroid search finished')
    return initialCentroids


def stationBelongTocluster(Config, CentroidList, StationMetaList):

    clusterList = []

    for i in StationMetaList:
        mind = 0
        c = 0

        for j in CentroidList:
            delta = abs(loc2degrees(j, i))
            c += 1
            if mind == 0:
                mind = delta
                i.member = c
            else:
                if delta < mind:
                    mind = delta
                    i.member = c
        clusterList.append(Station(i.net, i.sta, i.loc, i.comp, i.lat, i.lon,
                                   i.ele, i.dip, i.azi, i.gain, i.member))

    return clusterList


def calculateclusterCentre(Config, clusterStationList):

    newclusterVector = []
    cfg = ConfigObj(dict=Config)

    for i in range(1, cfg.Int('maxcluster')+1):
        sumlat = 0
        sumlon = 0
        clusterstationcounter = 0

        for j in clusterStationList:
            if i == j.member:
                sumlat += float(j.lat)
                sumlon += float(j.lon)
                clusterstationcounter += 1

        if clusterstationcounter == 0:
           newclusterVector.append(Centroid(0.0, 0.0, -1))
        else:
           scla = sumlat / clusterstationcounter
           sclo = sumlon / clusterstationcounter
           name = i
           newclusterVector.append(Centroid(scla,sclo,name))

    return newclusterVector

def compareclusterCentre(oldcluster, newcluster, Config):

    counter = 0
    for i in range(int(Config['maxcluster'])):
        if newcluster[i].rank == -1:
            continue

        delta = loc2degrees(oldcluster[i], newcluster[i])

        msg = str(i) + ' OLD: ' + str(oldcluster[i].lat) + ' ' + str(oldcluster[i].lon)
        msg+= ' NEW: ' + str(newcluster[i].lat) + ' ' + str(newcluster[i].lon) + ' DELTA: ' + str(delta)
        Logfile.add(msg)

        if delta < float(Config['centroidmindistance']):
           counter +=1
    return counter

def deleteFarStations(CentroidList, StationclusterList, Config):

    cfg = ConfigObj(dict=Config)
    stationdistance = int(cfg.Distance('stationdistance'))

    for i in CentroidList:
        for j in StationclusterList:
            if i.rank == j.member:
                if loc2degrees(i, j) > stationdistance:
                   j.member = -1

    for index,k in enumerate(StationclusterList):
        if k.member == -1:
            del StationclusterList[index]

    return StationclusterList


def filterclusterStationMinimumNumber(CentroidList, StationclusterList, Config):

    newCentroidList  = []
    newStationclusterList = []

    for i in CentroidList:
        counter = 0

        for j in StationclusterList:
            if i.rank == j.member:
                counter +=1
                streamID = j.net+'.'+j.sta+'.'+j.loc+'.'+j.comp

        if counter < int(Config['minclusterstation']): s1 = 'OUT'
        else:
           s1 = 'IN '
           newCentroidList.append(i)

        Logfile.red('Centroid %s %d %s %.2f %5.2f' % (i.rank, counter, s1,
                                                      i.lat, i.lon))

    for i in newCentroidList:
        for j in StationclusterList:
            if i.rank == j.member:
                newStationclusterList.append(j)

    return newStationclusterList, newCentroidList


def calcMeanCentroidDistance(CentroidList):

    sumdelta = 0

    for i in CentroidList:
        for j in CentroidList:
            sumdelta += loc2degrees(i, j)

    meanCentroidDistance = sumdelta / len(CentroidList)
    return meanCentroidDistance


def calcMinValue(CentroidList):

    mCD = calcMeanCentroidDistance(CentroidList)

    sumdelta = 0

    for i in CentroidList:
        for j in CentroidList:
            delta = loc2degrees(i, j)
            x=(delta - mCD)
            sumdelta += math.pow( x, 2 )

    minval = sumdelta / len(CentroidList)

    return minval


def write4Plot(Config, Origin, StationclusterList, CentroidList, Folder, flag):
    plotfolder = 'plot-'+str(flag)
    p = os.path.join(Folder['cluster'],plotfolder)
    if os.access(p,os.F_OK) == False:
        os.makedirs(p)

    fobjorigin = open(os.path.join(p,'event.orig'),'w')
    fobjorigin.write(Origin['lat']+','+Origin['lon'])
    fobjorigin.close()

    dst = os.path.join(p,'plot_cluster.py')
    import palantiri
    path = palantiri.__path__
    src = os.path.join(path[0], 'skeleton/plot_cluster.py')
    shutil.copy(src, dst)

    fobjcentroid = open(os.path.join(p,'event.centroid'),'w')

    for i in CentroidList:
        fobjcentroid.write(str(i.rank)+' '+str(i.lat)+' '+str(i.lon)+'\n')

    fobjcentroid.close()

    fobjstation = open(os.path.join(p,'event.stations'),'w')

    for i in StationclusterList:
        if i.loc == '--':  i.loc=''

        streamID = i.net+'.'+i.sta+'.'+i.loc+'.'+i.comp
        fobjstation.write(streamID+' '+i.lat+' '+i.lon+' '+str(i.member)+'\n')

    fobjstation.close()

    fobjcentroid = open(os.path.join(p,'event.statistic'),'w')
    mCD = calcMeanCentroidDistance(CentroidList)
    mv = calcMinValue(CentroidList)

    fobjcentroid.write(str(mCD)+' '+str(mv)+' '+str(len(CentroidList))+' '+str(len(StationclusterList))+'\n')
    fobjcentroid.close()


def km(Config, FilterMeta, Folder, Origin, flag):

    ic = DeterminedInitialCentroids(Config, FilterMeta)

    if len(ic) == 0:
        return

    counter = 0
    kmean(Config, ic, FilterMeta, counter, Folder, Origin, flag)


def endcheck(inputCentroid,FilterMeta,Config,Folder,Origin,flag):

    FM = deleteFarStations(inputCentroid,FilterMeta,Config)
    FMM,NCL = filterclusterStationMinimumNumber(inputCentroid,FM,Config)
    write4Plot(Config, Origin, FMM, NCL, Folder, flag)


def kmean(Config, inputCentroid, FilterMeta, counter, Folder, Origin, flag):

    counter += 1
    Logfile.add('COUNTER ' + str(counter) + ' CUTOFF ' + Config['cutoff'])

    cfg = ConfigObj(dict=Config)

    scl = stationBelongTocluster(Config, inputCentroid, FilterMeta)

    acounter = 1
    for a in inputCentroid:
        for i in scl:
            if acounter == i.member:
                delta = loc2degrees(i, a)
                if delta > cfg.Float('initialstationdistance'):
                    i.member = -1

        acounter+=1

    if counter == cfg.UInt('cutoff'):
        endcheck(inputCentroid, FilterMeta, Config, Folder, Origin, flag)
        sys.exit()

    nsc = calculateclusterCentre(Config, scl)
    t = compareclusterCentre(inputCentroid, nsc, Config)

    Logfile.add('ITERATIONSTEP: ---> ' + str(counter) + ' <-----------------------------')

    while t < cfg.UInt('maxcluster'):
        Logfile.add('number of arrays in KMEAN: ' + str(t))
        kmean(Config, nsc, FilterMeta, counter, Folder, Origin, flag)



    endcheck(inputCentroid, FilterMeta, Config, Folder, Origin, flag)
    sys.exit()


class clusterMain(MainObj):

    def __init__(self):

        parser = OptionParser(usage="%prog -f Eventpath ")
        parser.add_option("-f", "--evpath", type="string", dest="evpath",
                          help="evpath")

        (options, args) = parser.parse_args()
        self.eventpath = options.evpath

        Basic.checkExistsDir(self.eventpath, isAbort=True)
        Globals.setEventDir(self.eventpath)
        MainObj.__init__(self, self, 'ok', 'cluster_run.log',
                         'cluster.log')

    def init(self):
        return Globals.init()

    def process(self):

        t = time.time()
        C = config.Config(self.eventpath)

        Config = C.parseConfig('config')
        cfg = ConfigObj(dict=Config)
        Origin = C.parseConfig('origin')
        if cfg.pyrocko_download() is True:
            if cfg.quantity() == 'displacement':
                disp = True
            else:
                disp = False
            Meta = readpyrockostations(self.eventpath, disp, cfg)
        elif cfg.colesseo_input() is True:
            scenario = guts.load(filename=cfg.colosseo_scenario_yml())
            scenario_path = cfg.colosseo_scenario_yml()[:-12]
            Meta = readcolosseostations(scenario_path)
            events = scenario.get_events()
            ev = events[0]
            Origin['strike'] = str(ev.moment_tensor.strike1)
            Origin['rake'] = str(ev.moment_tensor.rake1)
            Origin['dip'] = str(ev.moment_tensor.dip1)
            Origin['lat'] = str(ev.lat)
            Origin['lon'] = str(ev.lon)
            Origin['depth'] = str(ev.depth/1000.)

        else:
            Meta = readMetaInfoFile(self.eventpath)
        Folder = createFolder(self.eventpath)

        FilterMeta = filterStations(Meta, Config, Origin)

        km(Config, FilterMeta, Folder, Origin, t)


        return True

    def finish(self):
        pass


def MainProc():

    mainObj = clusterMain()
    mainObj.run()


if __name__ == "__main__":

    MainProc()
