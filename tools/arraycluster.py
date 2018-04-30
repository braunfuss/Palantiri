import fnmatch
import os
from optparse import OptionParser
import logging
from obspy.taup.taup import locations2degrees
from obspy.core.util.geodetics import gps2DistAzimuth
import random
import time
import shutil
import sys
import config
from ConfigParser import SafeConfigParser

sys.setrecursionlimit(1500)

logger = logging.getLogger('CT')
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

parser = OptionParser(usage="%prog -f Eventpath ")
parser.add_option("-f", "--evpath", type="string", dest="evpath", help="evpath")

(options, args) = parser.parse_args()

#if options.evpath == None:
#        parser.error("non existing eventpath")

#evpath = options.evpath


class Station(object):
    def __init__(self, net, sta, loc, comp,lat=0,lon=0,ele=0,dip=0,azi=0,gain=0,member=-1):
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
        
class Centroid(object):
    def __init__(self,lat,lon,rank):
        self.lat = lat
        self.lon = lon
        self.rank = rank


def init():
    cDict = {}
    parser = SafeConfigParser()
    parser.read('global.conf')
    for section_name in parser.sections():
        for name, value in parser.items(section_name):
            cDict[name]=value
    return cDict


def readEvent(path):
    logger.info('\033[31m Parsing Origin \033[0m \n')
    cDict = {}
    parser = SafeConfigParser()
    parser.read('KYRGYZSTAN_2008-10-05_15-52-52.000000.origin')
    for section_name in parser.sections():
        for name, value in parser.items(section_name):
            cDict[name]=value
    return cDict

def readConfig(path):
    logger.info('\033[31m Parsing Config \033[0m \n')
    cDict = {}
    parser = SafeConfigParser()
    parser.read('KYRGYZSTAN_2008-10-05_15-52-52.000000.config')
    for section_name in parser.sections():
        for name, value in parser.items(section_name):
            cDict[name]=value
    return cDict

def readMetaInfoFile(EventPath):

    logger.info('\033[31m Parsing MetaInfoFile \033[0m \n')
    try:
        for i in os.listdir(EventPath):
            if fnmatch.fnmatch(i, '*.meta'):
                evfile = os.path.join(EventPath,i)

        MetaL = []
    
        fobj = open(evfile,'r')
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

            if fnmatch.fnmatch(comp,'*HZ'):
                MetaL.append(Station(net,sta,loc,comp,lat,lon,ele,dip,azi,gain))

        logger.info('\033[31m %d ENTRIES IN METAFILE FOUND \033[0m \n' % (len(MetaL)))
    except:
        logger.info('\033[31m METAFILE NOT READABLE \033[0m \n')


    return MetaL


def createFolder(EventPath):
    
    Folder = {}
    logger.info('\033[31m Create working environment \033[0m \n')
    if os.access(os.getcwd(), os.W_OK):
        
        basedir = os.path.join(EventPath,'work')
        sembdir = os.path.join(basedir, 'semblance')
        ascdir  = os.path.join(basedir, 'asc')
        mseeddir  = os.path.join(basedir, 'mseed')
        
        Folder['base'] = basedir
        Folder['semb'] = sembdir
        Folder['asc']  = ascdir
        Folder['mseed'] = mseeddir
        Folder['event'] = EventPath
        
    else:
        print"no write permissions"
        
    Folder['config'] = os.path.join(os.getcwd(),'..','skeleton')
    return Folder

def filterStations(StationList,Config,Origin):
    
    F = []
    minDist = int(Config['mindist'])
    maxDist = int(Config['maxdist'])
    
    o_lat = float(Origin['lat'])
    o_lon = float(Origin['lon'])
       
    logger.info('\033[31m Filter stations with configured parameters \033[0m')
    
    for i in StationList:

        sdelta = locations2degrees(o_lat, o_lon, float(i.lat), float(i.lon))

        if sdelta > minDist and sdelta < maxDist:
                    F.append(Station(i.net,i.sta,i.loc,i.comp,i.lat,i.lon,i.ele,i.dip,i.azi,i.gain))

    logger.info('\033[31m %d STATIONS LEFT IN LIST \033[0m'% len(F))
    return F

def checkStationAroundInitialCentroid(station,Config,StationMetaList):
    
    counter = 0
    for i in StationMetaList:
        sdelta = locations2degrees(float(station.lat), float(station.lon), float(i.lat), float(i.lon))
        if sdelta < int(Config['initialstationdistance']):
            counter +=1
    
    return counter

def addOK(station,stationList,Config,MetaList):
    
    t=0
    for i in stationList:
        sdelta = locations2degrees(float(station.lat), float(station.lon), float(i.lat), float(i.lon))
        
        if sdelta > float(Config['centroidmindistance']):
            aroundcounter = checkStationAroundInitialCentroid(station,Config,MetaList)
            if aroundcounter >= int(Config['minstationaroundinitialcluster']):
                t=1
            else:
                t=0
                return t
        else:
            t=0
            return t

    return t

def alreadyUsedIndex(index,usedIndexList):
    
    for i in usedIndexList:
        if i == index:
            return 1
        else:
            return 0


def createRandomInitialCentroids(Config,StationMetaList):
    
    print 'find initial cluster centers'
    
    initialCentroids = []
    usedIndexes = []
    random.seed(time.clock())
    
    counter = 1

    while len(initialCentroids) != int(Config['maxcluster']):

            randomIndex = random.randint(0, len(StationMetaList)-1)
            flag = alreadyUsedIndex(randomIndex,usedIndexes)
            usedIndexes.append(randomIndex)
        
            if len(usedIndexes) == len(StationMetaList):
                #sys.exit()
                pass
        
            counter += 1
            if len(initialCentroids) == 0:
                initialCentroids.append(StationMetaList[randomIndex])
            else:
                t = addOK(StationMetaList[randomIndex],initialCentroids,Config,StationMetaList)
                if t == 1:
                    initialCentroids.append(StationMetaList[randomIndex])
                    print 'found centroids til now ',len(initialCentroids)
                else:
                    continue

    return initialCentroids


def stationBelongToCluster(Config,CentroidList,StationMetaList):
    
    ClusterList = []
    
    for i in StationMetaList:
        mind = 100000
        c = 0
        for j in CentroidList:
            delta = locations2degrees(float(j.lat), float(j.lon), float(i.lat), float(i.lon))
            c+=1
            
            if mind > delta:
                mind=delta
                i.member = c
            
        ClusterList.append(Station(i.net,i.sta,i.loc,i.comp,i.lat,i.lon,i.ele,i.dip,i.azi,i.gain,i.member))
        
    
    return ClusterList


def calculateClusterCentre(Config,ClusterStationList):
    
    newClusterVector = []
    for i in range (1,int(Config['maxcluster'])+1):
        sumlat = 0
        sumlon = 0
        clusterstationcounter = 0
        for j in ClusterStationList:
            if i == j.member:
                sumlat += float(j.lat)
                sumlon += float(j.lon)
                clusterstationcounter += 1
                
        scla=sumlat/clusterstationcounter
        sclo=sumlon/clusterstationcounter
        name = i
        newClusterVector.append(Centroid(scla,sclo,name))
    
    return newClusterVector

def compareClusterCentre(oldCluster,newCluster,Config):
    
    counter = 0
    for i in range (int(Config['maxcluster'])):
        delta = locations2degrees(float(oldCluster[i].lat),float(oldCluster[i].lon), float(newCluster[i].lat),float(newCluster[i].lon))
        print i,' OLD: ',oldCluster[i].lat,oldCluster[i].lon,' NEW: ',newCluster[i].lat,newCluster[i].lon,' DELTA: ',delta
        if delta < float(Config['comparedelta']):
            print 'JO'
            counter +=1
        else:
            print 'NO'

    return counter

def deleteFarStations(CentroidList,StationClusterList,Config):
    
    for i in CentroidList:
        for j in StationClusterList:
            if i.rank == j.member:
                delta = locations2degrees(float(i.lat),float(i.lon), float(j.lat),float(j.lon))
                if delta > int(Config['stationdistance']):
                    j.member = -1
    
    #for index,k in enumerate(StationClusterList):
#        if k.member == -1:
 #           del StationClusterList[index]
            
    return StationClusterList


def filterClusterStationMinimumNumber(CentroidList,StationClusterList,Config):
    
    newCentroidList = []
    newStationClusterList = []
    
    for i in CentroidList:
        counter = 0
        for j in StationClusterList:
            if i.rank == j.member:
                counter+=1
                streamID = j.net+'.'+j.sta+'.'+j.loc+'.'+j.comp
                delta = locations2degrees(float(i.lat),float(i.lon), float(j.lat),float(j.lon))
                #print i.lat,i.lon,': ',streamID,j.lat,j.lon,delta
        
        if counter < int(Config['minclusterstation']):
            print 'Centroid ',i.rank, counter, ' OUT', i.lat,i.lon
        else:
            print 'Centroid ',i.rank, counter, ' IN', i.lat,i.lon
            newCentroidList.append(i)
    
#    for i in newCentroidList:
 #       for j in StationClusterList:
  #          if i.rank == j.member:
   #             newStationClusterList.append(j)

    #return newStationClusterList,newCentroidList
    return StationClusterList,newCentroidList


def write4Plot(Config,Origin,StationClusterList,CentroidList,Folder):
    
    p = os.path.join(Folder['event'],'plot')
    if os.access(p,os.F_OK) == False:
        os.makedirs(p)
    

    fobjorigin = open(os.path.join(p,'event.orig'),'w')
    fobjorigin.write(Origin['lat']+','+Origin['lon'])
    fobjorigin.close()
    
    fobjcentroid = open(os.path.join(p,'event.centroid'),'w')

    for i in CentroidList:
        fobjcentroid.write(str(i.rank)+' '+str(i.lat)+' '+str(i.lon)+'\n')

    fobjcentroid.close()
    
    fobjstation = open(os.path.join(p,'event.stations'),'w')
    for i in StationClusterList:
        if i.loc == '--':
            i.loc=''
        streamID = i.net+'.'+i.sta+'.'+i.loc+'.'+i.comp
        #print streamID,i.member
        fobjstation.write(streamID+' '+i.lat+' '+i.lon+' '+str(i.member)+'\n')
    fobjstation.close()
    
    
def km(Config,FilterMeta,Folder,Origin):
    
    
    ic = createRandomInitialCentroids(Config,FilterMeta)
    #sys.exit()
    
    counter = 0
    kmean(Config,ic,FilterMeta,counter,Folder,Origin)

def endcheck(inputCentroid,FilterMeta,Config,Folder,Origin):
    
    FM = deleteFarStations(inputCentroid,FilterMeta,Config)
    FMM,NCL = filterClusterStationMinimumNumber(inputCentroid,FM,Config)
    print len(FMM),'    ',len(FM)
    write4Plot(Config,Origin,FMM,NCL,Folder)

def kmean(Config,inputCentroid,FilterMeta,counter,Folder,Origin):
    
    counter += 1
    if counter == int(Config['cutoff']):
        
        #FM = deleteFarStations(inputCentroid,FilterMeta,Config)
        #FMM,NCL = filterClusterStationMinimumNumber(inputCentroid,FM,Config)
        #write4Plot(Config,Origin,FMM,inputCentroid,Folder)
        endcheck(inputCentroid,FilterMeta,Config,Folder,Origin)
        sys.exit()
        
        
    scl = stationBelongToCluster(Config,inputCentroid,FilterMeta)
    nsc = calculateClusterCentre(Config,scl)
    t = compareClusterCentre(inputCentroid,nsc,Config)
    
    print 'ITERATIONSTEP: ---> ',counter,' <----------------------------------------------\n\n'
    
    
    while t != int(Config['maxcluster']):
        print 'ANZAHL DER JO in KMEAN: ',t
        kmean(Config,nsc,FilterMeta,counter,Folder,Origin)
        
    #write4Plot(Config,Origin,FilterMeta,inputCentroid,Folder)
    endcheck(inputCentroid,FilterMeta,Config,Folder,Origin)
    sys.exit()
    

def run():
    
    eventpath = os.getcwd()
    Origin = readEvent(eventpath)
    Config = readConfig(eventpath)
    Meta   = readMetaInfoFile(eventpath)
    Folder = createFolder(eventpath)
    FilterMeta = filterStations(Meta,Config,Origin)
    
    
    km(Config,FilterMeta,Folder,Origin)
    

if __name__ == "__main__":
    
    run()