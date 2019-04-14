import os
import sys


sys.path.append('../tools/')
sys.path.append('../Common/')

from   optparse import OptionParser
import logging
import shutil
import fnmatch

import  Basic
import  Globals
import  Logfile
from ObspyFkt import loc2degrees, obs_kilometer2degrees

from config import Config
from cluster2 import Centroid,Station
from ConfigFile import ConfigObj

parser = OptionParser(usage="%prog -f eventpath ")
parser.add_option("-f", "--evpath", type="string", dest="evpath", help="evpath")

(options, args) = parser.parse_args()


class Result(object):

    def __init__(self, meanvalue, minvalue, centroidcount, usedstationcount,
                 path):

        self.meanvalue = meanvalue
        self.minvalue = minvalue
        self.centroidcount = centroidcount
        self.usedstationcount = usedstationcount
        self.path = path


class BestSolution(object):

    def __init__(self, station, cluster, lat, lon):

        self.station = station
        self.cluster = cluster
        self.lat = lat
        self.lon = lon


def getStatistics(clusterresults):
    p = os.path.join(options.evpath, 'stat.dat')
    fobjstat = open(p, 'w')
    resDict = {}
    for root, dirs, files in os.walk(clusterresults):
        for i in files:
            if i == 'event.statistic':
                fname = os.path.join(root, i)
                bspath = os.path.join('/', *fname.split('/')[:-1])
                fobj = open(fname, 'r')

                for line in fobj:
                    line = line.split()
                    resDict[fname] = Result(line[0], line[1], line[2], line[3],
                                            bspath)
                    fobjstat.write(('%s:  %s  %s  %s  %s\n') % (fname, line[0],
                                                                line[1],
                                                                line[2],
                                                                line[3]))

                fobj.close()

    fobjstat.close()
    return resDict


def getBestSolution(resultDictionary):

    bestsolution = -100000

    for i in resultDictionary.keys():
        resultDictionary[i].meanvalue = float(resultDictionary[i].meanvalue)
        if bestsolution < resultDictionary[i].meanvalue:
            bestsolution = resultDictionary[i].meanvalue

    L = []
    for j in resultDictionary.keys():
        if bestsolution == resultDictionary[j].meanvalue:
            L.append(resultDictionary[j])

    return L[0]


def copyCluster2EventConfig(ClusterDict, evpath):

    epath = os.path.join('/', *evpath.split('/')[:-1])
    t = epath.split(os.path.sep)[-1]
    fname = t+'.config'
    fullfname = os.path.join(epath, fname)
    L = []
    fobj = open(fullfname, 'r')

    for index, line in enumerate(fobj):
        L.append(line)

        if fnmatch.fnmatch(line, '*array parameter*'):
            firstend = index+1

        if fnmatch.fnmatch(line, '*beamforming method*'):
            secondbegin = index

    fobj.close()

    Confpart1 = L[:firstend]
    Confpart2 = L[secondbegin:]

    fobj = open(fullfname, 'w')
    fobj.write(''.join(Confpart1))
    nlist = ''

    for i in ClusterDict.keys():
        if len(ClusterDict[i]) > 0:
            nlist += 'r'+str(i)+','

    fobj.write(('networks=%s\n') % (nlist[:-1]))
    fobj.write('\n')

    for i in ClusterDict.keys():
        if len(ClusterDict[i]) > 0:
            aname = 'r'+str(i)
            fobj.write(('%s=%s\n') % (aname, ClusterDict[i][:-1]))
            fobj.write(('%srefstation=\n') %( aname))
            fobj.write(('%sphase=P\n') % (aname))

    fobj.write('\n')
    fobj.write(''.join(Confpart2))
    fobj.close()


def printBestSolution(solution):

    maxline = -100
    L = []
    M = []

    Logfile.add('eventpath: ', os.path.join(solution.path, 'event.stations'))
    fobj = open(os.path.join(solution.path, 'event.stations'), 'r')

    for line in fobj:
        line = line.split()
        M.append(int(line[3]))
        L.append(BestSolution(line[0], line[3], line[1], line[2]))

    fobj.close()

    maxline = max(M)

    C = {}
    fobjarrays = open(os.path.join(os.path.join('/',*solution.path.split('/')[:-2]),'arraycluster.dat'),'w')

    for i in range(int(maxline)+1):
        array = ''

        for j in L:
            if j.cluster == str(i):
                array+=j.station+'|'

        ar = 'r'+str(i)+'='+array[:-1]+'\n'
        C[i] = array
        fobjarrays.write(ar)

    fobjarrays.close()

    return C


def copyAndShowBestSolution(solution):

    src = os.path.join('/', *solution.path.split('/')[:-4])
    src = os.path.join(src, 'skeleton', 'clusterplot.sh')


def filterBestSolution(solution):

    evp = os.path.join('/', *solution.path.split('/')[:-2])
    C = Config(evp)
    Conf = C.parseConfig('config')
    cfg = ConfigObj(dict=Conf)

    SL = []
    M = []
    fobj = open(os.path.join(solution.path, 'event.stations'), 'r')

    for s in fobj:
        try:
            line = s.split()
            net, sta, loc, comp = line[0].split('.')

            slat = line[1]
            slon = line[2]
            smember = line[3]

            M.append(smember)
            SL.append(Station(net, sta, loc, comp, lat=slat, lon=slon,
                              member=smember))

        except Exception:
            Logfile.exception('filterBestSolution', '<' + s + '>')
            continue

    fobj.close()

    M = list(set(M))

    Logfile.add('number of clusters ' + str(len(M)),
                'number of stations ' + str(len(SL)))

    kd = obs_kilometer2degrees(cfg.Distance('intraclusterdistance'))
    Logfile.add('icd ' + str(kd))

    maxdist = -1

    for i in SL:
        counter = 0

        for k in SL:
            if i.member == '8' and k.member == '8':
                if i.getName() != k.getName():
                    delta = loc2degrees(i, k)

                    if delta > maxdist:
                        maxdist = delta

                    if delta < kd:
                        counter +=1


if __name__ == "__main__":
    rD = getStatistics(options.evpath)
    bs = getBestSolution(rD)
    CD = printBestSolution(bs)
    copyAndShowBestSolution(bs)

    copyCluster2EventConfig(CD, options.evpath)
