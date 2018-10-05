import os
import sys
import platform

# add local directories to import path

sys.path.append('../tools/')
sys.path.append('../Common/')

import cPickle as pickle

from config import Station

import fnmatch
import logging
import math
from   math  import sin, cos, atan2
import time
import subprocess

#    Import from common

import Basic
import Logfile
from   DataTypes  import Location
from   ConfigFile import ConfigObj
from   ObspyFkt   import loc2degrees, obs_TravelTimes
from pyrocko import cake
import numpy as np
km = 1000.
# -------------------------------------------------------------------------------------------------

logger = logging.getLogger(sys.argv[0])

d2r = math.pi/180.
r2d = 1./d2r

# -------------------------------------------------------------------------------------------------

class GridElem(object):
    def __init__(self, lat, lon, depth, tt, delta):
        self.lat   = lat
        self.lon   = lon
        self.depth = depth
        self.tt= tt
        self.delta = delta

# -------------------------------------------------------------------------------------------------

class TTTGrid(object):
    def __init__(self, dimZ, mint,maxt,Latul,Lonul,Lator,Lonor,GridArray):
        self.dimZ  = dimZ
        self.mint  = mint
        self.maxt  = maxt
        self.Latul = Latul
        self.Lonul = Lonul
        self.Lator = Lator
        self.Lonor = Lonor
        self.GridArray = GridArray

# -------------------------------------------------------------------------------------------------

class MinTMaxT(object):
    def __init__(self, mint,maxt):
        self.mint = mint
        self.maxt = maxt

# -------------------------------------------------------------------------------------------------

def filterStations(StationList,Config,Origin,network):

    F = []
    cfg = ConfigObj(dict=Config)

    minDist, maxDist = cfg.FloatRange('mindist', 'maxdist')
    origin = Location(Origin['lat'], Origin['lon'])

    Logfile.red('Filter stations with configured parameters...')
    print 'nr networks = ', len(network)
    print 'nr stations = ', len(StationList)

    for j in network:
        for i in StationList:
            if fnmatch.fnmatch(i.getcmpName(),j) :
                pos= Location(i.lat, i.lon)
                sdelta = loc2degrees(origin, pos)
                if sdelta > minDist and sdelta < maxDist:
                    s = Station(i.net,i.sta,i.loc,i.comp,i.lat,i.lon,i.ele,i.dip,i.azi,i.gain)
                    if s not in F: F.append(s)
    Logfile.red('%d STATIONS LEFT IN LIST' % len(F))
    return F

# -------------------------------------------------------------------------------------------------

def calctakeoff(Station,Event,Config):

    de   = loc2degrees(Event, Station)
    Phase = cake.PhaseDef('P')
    model = cake.load_model()
    arrivals= model.arrivals([de,de], phases=Phase, zstart=Event.depth*km)

    return arrivals[0].takeoff_angle()

# -------------------------------------------------------------------------------------------------

def bearing(Station ,Event):

        #Convert to radians.
        lat1 = d2r * float(Station.lat)
        lon1 = d2r * float(Station.lon)
        lat2 = d2r * float(Event.lat)
        lon2 = d2r * float(Event.lon)

        #Compute the angle.
        x = sin(lon1-lon2) * cos(lat2);
        y = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1-lon2);

        angle = -atan2(x,y);

        if(angle < 0.0 ) :  angle  += math.pi * 2.0;

        angle = r2d * angle
        return angle;
# -------------------------------------------------------------------------------------------------

def backazi(Station, Event):

        lat1 = d2r * float(Station.lat)
        lon1 = d2r * float(Station.lon)
        lat2 = d2r * float(Event.lat)
        lon2 = d2r * float(Event.lon)

        #Compute the angle.
        x = sin(lon1-lon2) * cos(lat2)
        y = cos(lat1)      * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1-lon2)

        angle = -atan2(x,y);

        if(angle < 0.0 ):
            angle  += math.pi * 2.0;

        #And convert result to degrees.
        angle = r2d * angle

        return angle

# -------------------------------------------------------------------------------------------------
def calcTTTAdv(Config,station,Origin,flag,arrayname,Xcorrshift=None,Refshift=None):

    phasename =('%sphase') %(os.path.basename(arrayname))

    cfg= ConfigObj(dict=Config)
    dimX= cfg.Int('dimx')
    dimY= cfg.Int('dimy')
    gridspacing = cfg.Float('gridspacing')

    o_lat   = float(Origin['lat'])
    o_lon   = float(Origin['lon'])
    o_depth = float(Origin['depth'])

    oLator = o_lat + dimX/2;   oLonor = o_lon + dimY/2
    oLatul = 0
    oLonul = 0
    o_dip = 80.
    plane = False

    TTTGridMap = {}
    LMINMAX= []
    GridArray  = {}
    locStation = Location (station.lat, station.lon)
    sdelta = loc2degrees(Location(o_lat, o_lon), locStation)
    Phase = cake.PhaseDef(Config[phasename])
    model = cake.load_model()

    z = 0
    if plane is True:
        depth = np.linspace(0., 40., num=dimY)
        for i in xrange(70):
            oLatul = o_lat-((dimX-1)/2) * gridspacing + i * gridspacing
            if z == 0 and i == 0:
                Latul = oLatul
            o = 0
            start_time = time.clock()

            for j in xrange(40):
                oLonul = o_lon -((dimY-1)/2)* gridspacing + j * gridspacing/np.cos(o_dip)
                if o==0 and j==0: Lonul = oLonul
                de = loc2degrees(Location(oLatul, oLonul), locStation)
                arrivals = model.arrivals([de,de], phases=Phase, zstart=depth[j]*km, zstop=0.)
            try:
                ttime = arrivals[0].t
            except:
                try:
                    arrivals = model.arrivals([de,de], phases=Phase, zstart=depth[j]*km-2.5, zstop=depth[j]*km+2.5, refine=True)
                    ttime = arrivals[0].t
                except:
                    tt = obs_TravelTimes(de, o_depth)
                    for k in tt:
                        if k['phase_name'] == 'P' or k['phase_name'] ==('%sdiff')%(Config[phasename]):
                            ttime = k ['time']
                        print "Something wrong with phase arrival too large\
                             distances choosen?"

                GridArray[(i,j)] = GridElem(oLatul, oLonul, depth[j],ttime,de)
                LMINMAX.append(ttime)
                if int(Config['xcorr']) == 1:
                    ttime = ttime-float(Xcorrshift[station.getName()].shift)\
                     - Refshift
                GridArray[(i, j)] = GridElem(oLatul, oLonul, o_depth, ttime, de)
                LMINMAX.append(ttime)

                if ttime == 0:
                    print '\033[31mAvailable phases for station %s in range %f deegree\033[0m'%(station,de)
                    print '\033[31m'+'|'.join([str(item['phase_name']) for item in tt])+'\033[0m'
                    print '\033[31myou tried phase %s\033[0m'%(Config[phasename])
                    raise Exception("\033[31mILLEGAL: phase definition\033[0m")
    else:

        for i in xrange(dimX):
            oLatul = o_lat -((dimX-1)/2) * gridspacing + i * gridspacing

            if z == 0 and i == 0 :
                Latul = oLatul
            o=0
            for j in xrange(dimY):
                oLonul = o_lon-((dimY-1)/2) * gridspacing + j * gridspacing

                if o==0 and j==0: Lonul = oLonul
                de = loc2degrees(Location(oLatul, oLonul), locStation)
                arrivals = model.arrivals([de, de], phases=Phase,
                                          zstart=o_depth*km)
            try:
                ttime = arrivals[0].t
            except:
                try:
                    arrivals = model.arrivals([de, de], phases=Phase,
                                              zstart=o_depth*km-2.5,
                                              zstop=o_depth*km+2.5,
                                              refine=True)
                    ttime = arrivals[0].t
                except:
                    tt = obs_TravelTimes(de, o_depth)
                    l = tt[0]
                    ttime = l.time

                GridArray[(i, j)] = GridElem(oLatul, oLonul, o_depth,ttime,de)
                LMINMAX.append(ttime)
                if int(Config['xcorr']) == 1:
                    ttime = ttime-float(Xcorrshift[station.getName()].shift)\
                     - Refshift
                GridArray[(i, j)] = GridElem(oLatul, oLonul, o_depth, ttime, de)
                LMINMAX.append(ttime)

                if ttime == 0:
                    print '\033[31mAvailable phases for station %s in range %f deegree\033[0m'%(station,de)
                    print '\033[31m'+'|'.join([str(item['phase_name']) for item in tt])+'\033[0m'
                    print '\033[31myou tried phase %s\033[0m'%(Config[phasename])

                    raise Exception("\033[31mILLEGAL: phase definition\033[0m")

    mint = min(LMINMAX)
    maxt = max(LMINMAX)
    TTTGridMap[station.getName()] = TTTGrid(o_depth,mint,maxt,Latul,Lonul,oLator,oLonor,GridArray)
    k = MinTMaxT(mint,maxt)

    Basic.dumpToFile(str(flag)  + '-ttt.pkl', TTTGridMap)
    Basic.dumpToFile('minmax-'  + str(flag) + '.pkl', k)
    Basic.dumpToFile('station-' + str(flag) + '.pkl', station)

# -------------------------------------------------------------------------------------------------

def calcTTTAdvTauP(Config,station,Origin,flag,Xcorrshift=None,Refshift=None):

    cfg = ConfigObj(dict=Config)

    dimX= cfg.Int('dimx')
    dimY= cfg.Int('dimy')
    gridspacing = cfg.Float('gridspacing')

    o_lat   = float(Origin['lat'])
    o_lon   = float(Origin['lon'])
    o_depth = float(Origin['depth'])

    oLator = o_lat + dimX/2
    oLonor = o_lon + dimY/2
    oLatul = 0
    oLonul = 0

    TTTGridMap = {}
    LMINMAX= []
    GridArray  = {}
    locStation = Location(station.lat, station.lon)

    sdelta = loc2degrees(Location(o_lat, o_lon), locStation)
    Logfile.add('TTT PROCESS %d STATION: %s --> DELTA: %f'%(flag,station.getName(),sdelta))

    inputpath  = str(flag)+'-'+station.getName()+".input";
    outputpath = str(flag)+'-'+station.getName()+".output";
    errorpath  = str(flag)+'-'+station.getName()+'.error'

    fobjinput = open(inputpath,'w')

    fobjinput.write('s\n')
    fobjinput.write(('%s %s\n')%(station.lat,station.lon))
    fobjinput.write('h\n')
    fobjinput.write(('%s\n')%(o_depth))

    for i in xrange(dimX):
        oLatul = o_lat -((dimX-1)/2) * gridspacing + i * gridspacing

        for j in xrange(dimY):
            oLonul = o_lon -((dimY-1)/2) * gridspacing + j * gridspacing

            fobjinput.write('e\n')
            fobjinput.write(('%s %s\n')%(oLatul,oLonul))
    # endfor

    fobjinput.close()

    cmd =('taup_time -ph P -mod ak135 -time -o %s < %s > %s') %(outputpath,inputpath,errorpath)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()

    L  = []
    output = open(outputpath,'r')
    'OUTPUT: ', outputpath

    for k in output:
        k = k.split()

        if len(k) == 1:
            tt = k[0].replace('\n','')
            tt = float(tt)-float(Xcorrshift[station.getName()].shift)
            L.append(tt)
    # endfor

    output.close()

    z=0

    for i in xrange(dimX):
            oLatul = o_lat -((dimX-1)/2) * gridspacing + i * gridspacing

            if z == 0 and i == 0:
                Latul = oLatul
            o=0

            for j in xrange(dimY):
                oLonul = o_lon -((dimY-1)/2) * gridspacing + j * gridspacing

                if o==0 and j==0:
                    Lonul = oLonul

                de = loc2degrees(Location(oLatul, oLonul), locStation)
                time = L[i * dimX+j]

                GridArray[(i,j)] = GridElem(oLatul, oLonul, o_depth,time,de)
                LMINMAX.append(time)
    # endfor

    mint = float(min(LMINMAX))
    maxt = float(max(LMINMAX))
    k = MinTMaxT(mint,maxt)

    TTTGridMap[station.getName()] = TTTGrid(o_depth,mint,maxt,Latul,Lonul,oLator,oLonor,GridArray)

    tttname = str(flag)+'-ttt.pkl'
    Basic.dumpToFile(tttname, TTTGridMap)
    Basic.dumpToFile('minmax-'+str(flag)+'.pkl', k)

    try:
        os.remove(inputpath)
        os.remove(outputpath)
        os.remove(errorpath)

    except:
        Logfile.exception('cannot delete files')

# -------------------------------------------------------------------------------------------------

def calcak135parameter(Event):

    vp = 0
    vs = 0
    rho= 0
    tmprho = 0
    tmpvp  = 0
    tmpvs  = 0
    depth  = "%.02f" % float(Event.depth)
    Logfile.add('MODELDEPTH ' + depth)

    Basic.checkFileExists('ak135.model', isAbort=True)

    try:
        fobj = open('ak135.model','r')
        index = {}

        for counter,line in enumerate(fobj):
            line = line.split()
            index[counter] = line[0]

        L = []

        for i in index.iterkeys():
            if fnmatch.fnmatch(index[i], str(depth)+'*'):
                L.append(i)

        fobj.seek(0)

        if len(L) == 0:
            for counter,line in enumerate(fobj):
                line = line.split()

                if float(depth) > float(line[0]):
                    x = counter

            fobj.seek(0)
            begin = x
            end = x+1
            I = []

            for counter, line in enumerate(fobj):
                if counter == begin:  I.append(line)
                if counter == end:    I.append(line)

            dist = float(I[1].split()[0]) - float(I[0].split()[0])
            rhodist = float(I[1].split()[1]) - float(I[0].split()[1])
            vpdist = float(I[1].split()[2]) - float(I[0].split()[2])
            vsdist = float(I[1].split()[3]) - float(I[0].split()[3])
            depthdist = float(depth) - float(I[0].split()[0])

            rhofact = rhodist/dist
            vpfact = vpdist/dist
            vsfact = vsdist/dist

            rho = float(I[0].split()[1]) + rhofact * depthdist
            vp = float(I[0].split()[2]) + vpfact * depthdist
            vs = float(I[0].split()[3]) + vsfact * depthdist

        # endif

        n = len(L)

        if n == 1:
            for counter, line in enumerate(fobj):
                for a in L:
                    if counter == a:
                        line = line.split()
                        tmprho = float(line[1])
                        tmpvp = float(line[2])
                        tmpvs = float(line[3])
            # endfor

            rho = tmprho / n
            vp = tmpvp  / n
            vs = tmpvs  / n
        # endif

        if n == 2:
            for counter,line in enumerate(fobj):
                for a in L:
                    if counter == a:
                        line = line.split()
                        tmprho += float(line[1])
                        tmpvp += float(line[2])
                        tmpvs += float(line[3])
            # endfor

            rho = tmprho / n
            vp = tmpvp  / n
            vs = tmpvs  / n
        # endif

        fobj.close()

    except:
        Logfile.exception('Model File not found Exception')

    return vp, vs, rho

# -------------------------------------------------------------------------------------------------

def dubcup(rho, vp, vs, stri, dip, rak, azi, phi):

    dstri = float(stri) * d2r;
    ddip  = float(dip)  * d2r;
    drak  = float(rak)  * d2r;
    dazi  = float(azi)  * d2r;

    rad1  =  cos(drak) * sin(ddip) * sin(2.0 *(dazi-dstri));
    rad2  = -cos(drak) * cos(ddip) * cos(dazi-dstri);
    rad3  =  sin(drak) * sin(2.0 * ddip);

    rad4   =       -sin(dazi-dstri) * sin(dazi-dstri);
    rad5   =        sin(drak)       * cos(2.0 * ddip) * sin(dazi-dstri);
    rad6   =  0.5 * cos(drak)       * sin(ddip)       * sin(2.0 *(dazi-dstri));
    rad7   = -0.5 * sin(drak)       * sin(2.0 * ddip) *(1.0-rad4);
    ph = float(phi * d2r);
    radra1 = sin(2.0 * ph);

    #/* SV waves at source */

    radra2 = cos(2.0 * ph);
    ducusw = rad5 * radra2 + rad2 * radra2 + rad6 * radra1 + rad7 * radra1;

    #/* P waves at source */

    radra3 = sin(ph) * sin(ph);
    ducupw = rad1 * radra3 + rad2 * radra1 + rad3 *(cos(ph) * cos(ph) + radra3 * rad4) + rad5 * radra1;

    #/* SH waves at source */
    rad8   = cos(drak) * cos(ddip) * sin(dazi-dstri);
    rad9   = cos(drak) * sin(ddip) * cos(2.0 *(dazi-dstri));
    rad10  = sin(drak) * cos(2.0 * ddip) * cos(dazi-dstri);
    rad11  = -0.5 * sin(drak) * sin(2.0 * ddip) * sin(2.0 *(dazi-dstri));
    ducush = rad8 * cos(ph)   + rad9 * sin(ph)  + rad10 * cos(ph) + rad11 * sin(ph);

    if ducusw < 0.0:  svsign = -1.0;
    else:             svsign =  1.0;

    if ducupw < 0.0:  psign = -1.0;
    else:             psign =  1.0;

    if ducush < 0.0:  shsign = -1.0
    else:             shsign =  1.0

    pamp = ducupw;
    svamp = ducusw;
    shamp = ducush;

    return psign
'''
/* *pamp = ducupw*(4.0*PI*rho*vp*vp*vp);
 * *svamp = ducusw*(4.0*PI*rho*vs*vs*vs);
 * *shamp = ducush*(4.0*PI*rho*vs*vs*vs);
 */
 '''
