import os
import sys
import math
from math import radians, cos, sin, atan2
import logging
from obspy.core.utcdatetime import UTCDateTime
from palantiri.common import Basic, Globals, Logfile, DataTypes
from palantiri.common.DataTypes import Location
from palantiri.common.ObspyFkt import loc2degrees
from palantiri.common.ConfigFile import ConfigObj, OriginCfg, SynthCfg
import time
import numpy as num
from collections import OrderedDict, defaultdict
from pyrocko.gf import ws, LocalEngine, Target, DCSource, RectangularSource
from pyrocko import util, pile, model, catalog, gf, cake, io, trace
from pyrocko.guts import Object, String, Float, List
from palantiri.process.semp import otest
from palantiri.process.beam_stack import BeamForming
from pyrocko.gf import STF
from palantiri.process.stacking import PWS_stack
from palantiri.process import sembCalc
from palantiri.process import ttt, waveform, sembCalc, trigger
from palantiri.tools import config
from palantiri.process.array_crosscorrelation_v4  import Xcorr, cmpFilterMetavsXCORR, getArrayShiftValue
import cPickle  as pickle
import times

km = 1000.


logger = logging.getLogger('ARRAY-MP')


class CombiSource(gf.Source):
    '''Composite source model.'''

    discretized_source_class = gf.DiscretizedMTSource

    subsources = List.T(gf.Source.T())

    def __init__(self, subsources=[], **kwargs):
        if subsources:

            lats = num.array(
                [subsource.lat for subsource in subsources], dtype=num.float)
            lons = num.array(
                [subsource.lon for subsource in subsources], dtype=num.float)

            assert num.all(lats == lats[0]) and num.all(lons == lons[0])
            lat, lon = lats[0], lons[0]

            # if not same use:
            # lat, lon = center_latlon(subsources)

            depth = float(num.mean([p.depth for p in subsources]))
            t = float(num.mean([p.time for p in subsources]))
            kwargs.update(time=t, lat=float(lat), lon=float(lon), depth=depth)

        gf.Source.__init__(self, subsources=subsources, **kwargs)

    def get_factor(self):
        return 1.0

    def discretize_basesource(self, store, target=None):

        dsources = []
        t0 = self.subsources[0].time
        for sf in self.subsources:
            assert t0 == sf.time
            ds = sf.discretize_basesource(store, target)
            ds.m6s *= sf.get_factor()
            dsources.append(ds)

        return gf.DiscretizedMTSource.combine(dsources)




class CakeTiming(Object):
    '''Calculates and caches phase arrivals.
    :param fallback_time: returned, when no phase arrival was found for the
                        given depth-distance-phase-selection-combination

    E.g.:
    definition = 'first(p,P)-20'
    CakeTiming(definition)'''
    phase_selection = String.T()
    fallback_time = Float.T(optional=True)

    def __init__(self, phase_selection, fallback_time=None):
        self.arrivals = defaultdict(dict)
        self.fallback_time = fallback_time
        self.which = None
        self.phase_selection = phase_selection
        _phase_selection = phase_selection
        if '+' in _phase_selection:
            _phase_selection, self.offset = _phase_selection.split('+')
            self.offset = float(self.offset)
        elif '-' in _phase_selection:
            _phase_selection, self.offset = _phase_selection.split('-')
            self.offset = float(self.offset)
            self.offset = -self.offset

        if 'first' in _phase_selection:
            self.which = 'first'
        if 'last' in _phase_selection:
            self.which = 'last'
        if self.which:
            _phase_selection = self.strip(_phase_selection)

        self.phases = _phase_selection.split('|')

    def return_time(self, ray):
        if ray is None:
            return self.fallback_time
        else:
            return ray.t + self.offset

    def t(self, mod, z_dist, get_ray=False):
        ''':param phase_selection: phase names speparated by vertical bars
        :param z_dist: tuple with (depth, distance)
        '''
        z, dist = z_dist
        if (dist, z) in self.arrivals.keys():
            return self.return_time(self.arrivals[(dist, z)])

        phases = [cake.PhaseDef(pid) for pid in self.phases]
        arrivals = mod.arrivals(
            distances=[dist*cake.m2d], phases=phases, zstart=z)
        if arrivals == []:
            logger.warn(
                'no phase at d=%s, z=%s. (return fallback time)' % (dist, z))
            want = None
        else:
            want = self.phase_selector(arrivals)
        self.arrivals[(dist, z)] = want
        if get_ray:
            return want
        else:
            return self.return_time(want)

    def phase_selector(self, _list):
        if self.which == 'first':
            return min(_list, key=lambda x: x.t)
        if self.which == 'last':
            return max(_list, key=lambda x: x.t)

    def strip(self, ps):
        ps = ps.replace(self.which, '')
        ps = ps.rstrip(')')
        ps = ps.lstrip('(')
        return ps


class Timings(Object):
    timings = List.T(CakeTiming.T())

    def __init__(self, timings):
        self.timings = timings


class SembMax (object):
    '''
    class to store sembmax object for each grid point
    '''

    def __init__(self, lat, lon, semb):

        self.lat  = lat
        self.lon  = lon
        self.semb = semb

# -------------------------------------------------------------------------------------------------

class FileSembMax(object):
    '''
    class to strore sembmax object for the sembmaxvalue file
    '''
    def __init__(self,istep,sembmaxX,sembmaxY,sembmax,usedarrays,delta,azi,deltakm):

        self.istep  = istep
        self.sembmaxX   = sembmaxX
        self.sembmaxY   = sembmaxY
        self.sembmax= sembmax
        self.usedarrays = usedarrays
        self.delta  = delta
        self.azi= azi
        self.deltakm= deltakm

    def get(self):
        return ('%d %.2f %.2f %f %d %03f %f %03f\n' % (self.istep,self.sembmaxX,self.sembmaxY,
                                                       self.sembmax,self.usedarrays,self.delta,
                                                       self.azi,self.delta*119.19))

# -------------------------------------------------------------------------------------------------

def toAzimuth (latevent,lonevent,latsource,lonsource):
        '''
        method to calculate azimuth between two points
        '''
        # Convert to radians.
        lat1 = radians (latsource);
        lon1 = radians (lonsource);
        lat2 = radians (latevent);
        lon2 = radians (lonevent);

       # Compute the angle.
        x =  sin(lon1-lon2 ) * cos(lat2);
        y =  cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1-lon2);
        angle = -atan2 (x,y);

        if angle < 0.0 :
         angle  += math.pi * 2.0;

       #And convert result to degrees.
        angle = math.degrees (angle)
        angle = '%02f'%angle

        return angle;

# -------------------------------------------------------------------------------------------------

def writeSembMatricesSingleArray (SembList,Config,Origin,arrayfolder,ntimes,switch):
    '''
    method to write semblance matrizes from one processes to file for each timestep
    '''
    logger.info ('start write semblance matrices')

    cfg= ConfigObj (dict=Config)
    origin = OriginCfg (Origin)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    latv   = []
    lonv   = []

    gridspacing = cfg.Float ('gridspacing')
    migpoints   = dimX * dimY

    o_lat   = origin.lat()         # float (Origin['lat'])
    o_lon   = origin.lon()         # float (Origin['lon'])
    oLatul  = 0
    oLonul  = 0

    z=0

    for i in xrange(dimX):
         oLatul = o_lat - ((dimX-1)/2) * gridspacing + i*gridspacing

         if z == 0 and i == 0:
             Latul = oLatul
         o=0

         for j in xrange (dimY):
               oLonul = o_lon - ((dimY-1)/2) * gridspacing + j * gridspacing

               if o==0 and j==0:  Lonul = oLonul

               latv.append (oLatul)
               lonv.append (oLonul)
    #endfor

    rc  = UTCDateTime (Origin['time'])
    rcs = '%s-%s-%s_%02d:%02d:%02d'% (rc.day,rc.month,rc.year, rc.hour,rc.minute,rc.second)
    d   = rc.timestamp

    for a, i in enumerate(SembList):
        #logger.info('timestep %d' % a)

        fobj = open (os.path.join (arrayfolder,'%s-%s_%03d.ASC' % (switch,Origin['depth'],a)),'w')
        fobj.write ('# %s , %s\n' % (d,rcs))
        fobj.write ('# step %ds| ntimes %d| winlen: %ds\n' % (step,ntimes,winlen))
        fobj.write ('# \n')
        fobj.write ('# southwestlat: %.2f dlat: %f nlat: %f \n'%(Latul,gridspacing,dimX))
        fobj.write ('# southwestlon: %.2f dlon: %f nlon: %f \n'%(Lonul,gridspacing,dimY))
        fobj.write ('# ddepth: 0 ndepth: 1 \n')

        for j in range (migpoints):
            x= latv[j]
            y= lonv[j]
            z= origin.depth()         # float(Origin['depth'])
            semb = i[j]

            fobj.write ('%.2f %.2f %.2f %.20f\n' % (x,y,z,semb))
        #endfor

        fobj.close()
    #endfor

# -------------------------------------------------------------------------------------------------

def collectSemb (SembList,Config,Origin,Folder,ntimes,arrays,switch):
    '''
    method to collect semblance matrizes from all processes and write them to file for each timestep
    '''
    Logfile.add ('start collect in collectSemb')

    cfg= ConfigObj (dict=Config)
    origin = ConfigObj (dict=Origin)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    latv= []
    lonv= []

    gridspacing = cfg.Float ('gridspacing')
    migpoints   = dimX * dimY
    o_lat   = origin.lat()         # float (Origin['lat'])
    o_lon   = origin.lon()         # float (Origin['lon'])
    oLatul  = 0
    oLonul  = 0

    z=0

    for i in xrange(dimX):
         oLatul = o_lat - ((dimX-1)/2) * gridspacing + i*gridspacing

         if z == 0 and i == 0 :
             Latul = oLatul
         o=0

         for j in xrange (dimY):
               oLonul = o_lon - ((dimY-1)/2) * gridspacing + j*gridspacing

               if o==0 and j==0:
                    Lonul = oLonul

               latv.append (oLatul)
               lonv.append (oLonul)


    tmp=1
    for a in SembList:
        tmp *= a
    #sys.exit()

    sembmaxvaluev = num.ndarray (ntimes,dtype=float)
    sembmaxlatv   = num.ndarray (ntimes,dtype=float)
    sembmaxlonv   = num.ndarray (ntimes,dtype=float)

    rc= UTCDateTime(Origin['time'])
    rcs= '%s-%s-%s_%02d:%02d:%02d'% (rc.day,rc.month,rc.year, rc.hour,rc.minute,rc.second)
    d = rc.timestamp
    usedarrays = 5

    folder  = Folder['semb']
    fobjsembmax = open (os.path.join (folder,'sembmax_%s.txt' % (switch)),'w')
    for a, i in enumerate(tmp):
        logger.info('timestep %d' % a)


        fobj  = open (os.path.join (folder,'%s-%s_%03d.ASC' % (switch,Origin['depth'],a)),'w')
        #fobj = open (os.path.join (folder, '%03d.ASC'    % a),'w')

        fobj.write ('# %s , %s\n' % (d,rcs))
        fobj.write ('# step %ds| ntimes %d| winlen: %ds\n' % (step,ntimes,winlen))
        fobj.write ('# \n')
        fobj.write ('# southwestlat: %.2f dlat: %f nlat: %f \n'%(Latul,gridspacing,dimX))
        fobj.write ('# southwestlon: %.2f dlon: %f nlon: %f \n'%(Lonul,gridspacing,dimY))
        fobj.write ('# ddepth: 0 ndepth: 1 \n')


        sembmax  = 0
        sembmaxX = 0
        sembmaxY = 0

        origin = DataTypes.dictToLocation (Origin)
        uncert = num.std(i) #maybe not std?
        for j in range(migpoints):
            x= latv[j]
            y= lonv[j]
            semb = i[j]

            fobj.write ('%.2f %.2f %.20f\n' % (x,y,semb))

            if  semb > sembmax:
                sembmax  = semb;# search for maximum and position of maximum on semblance grid for given time step
                sembmaxX = x;
                sembmaxY = y;

        delta = loc2degrees (Location (sembmaxX, sembmaxY), origin)
        azi   = toAzimuth (float(Origin['lat']), float(Origin['lon']),float(sembmaxX), float(sembmaxY))

        sembmaxvaluev[a] = sembmax
        sembmaxlatv[a]   = sembmaxX
        sembmaxlonv[a]   = sembmaxY

        fobjsembmax.write ('%d %.2f %.2f %.20f %.20f %d %03f %f %03f\n' % (a*step,sembmaxX,sembmaxY,sembmax,uncert,usedarrays,delta,float(azi),delta*119.19))
        fobj.close()


    fobjsembmax.close()

    durationpath  = os.path.join (folder, "duration.txt")
    trigger.writeSembMaxValue (sembmaxvaluev,sembmaxlatv,sembmaxlonv,ntimes,Config,Folder)
    trigger.semblancestalta (sembmaxvaluev,sembmaxlatv,sembmaxlonv)

def collectSembweighted(SembList,Config,Origin,Folder,ntimes,arrays,switch, weights):
    '''
    method to collect semblance matrizes from all processes and write them to file for each timestep
    '''
    Logfile.add ('start collect in collectSemb')

    cfg= ConfigObj (dict=Config)
    origin = ConfigObj (dict=Origin)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    latv= []
    lonv= []

    gridspacing = cfg.Float ('gridspacing')
    migpoints   = dimX * dimY
    o_lat   = origin.lat()         # float (Origin['lat'])
    o_lon   = origin.lon()         # float (Origin['lon'])
    oLatul  = 0
    oLonul  = 0

    z=0

    for i in xrange(dimX):
         oLatul = o_lat - ((dimX-1)/2) * gridspacing + i*gridspacing

         if z == 0 and i == 0 :
             Latul = oLatul
         o=0

         for j in xrange (dimY):
               oLonul = o_lon - ((dimY-1)/2) * gridspacing + j*gridspacing

               if o==0 and j==0:
                    Lonul = oLonul

               latv.append (oLatul)
               lonv.append (oLonul)


    tmp=1
    for a, w in zip(SembList, weights):
        tmp *= a
    #sys.exit()

    sembmaxvaluev = num.ndarray (ntimes,dtype=float)
    sembmaxlatv   = num.ndarray (ntimes,dtype=float)
    sembmaxlonv   = num.ndarray (ntimes,dtype=float)

    rc= UTCDateTime(Origin['time'])
    rcs= '%s-%s-%s_%02d:%02d:%02d'% (rc.day,rc.month,rc.year, rc.hour,rc.minute,rc.second)
    d = rc.timestamp
    usedarrays = 5

    folder  = Folder['semb']
    fobjsembmax = open (os.path.join (folder,'sembmax_%s.txt' % (switch)),'w')

    for a, i in enumerate(tmp):
        logger.info('timestep %d' % a)


        fobj  = open (os.path.join (folder,'%s-%s_%03d._weighted_semblance.ASC' % (switch,Origin['depth'],a)),'w')
        #fobj = open (os.path.join (folder, '%03d.ASC'    % a),'w')

        fobj.write ('# %s , %s\n' % (d,rcs))
        fobj.write ('# step %ds| ntimes %d| winlen: %ds\n' % (step,ntimes,winlen))
        fobj.write ('# \n')
        fobj.write ('# southwestlat: %.2f dlat: %f nlat: %f \n'%(Latul,gridspacing,dimX))
        fobj.write ('# southwestlon: %.2f dlon: %f nlon: %f \n'%(Lonul,gridspacing,dimY))
        fobj.write ('# ddepth: 0 ndepth: 1 \n')


        sembmax  = 0
        sembmaxX = 0
        sembmaxY = 0

        origin = DataTypes.dictToLocation (Origin)
        uncert = num.std(i) #maybe not std?
        for j in range(migpoints):
            x= latv[j]
            y= lonv[j]
            semb = i[j]

            fobj.write ('%.2f %.2f %.20f\n' % (x,y,semb))

            if  semb > sembmax:
                sembmax  = semb;# search for maximum and position of maximum on semblance grid for given time step
                sembmaxX = x;
                sembmaxY = y;

        delta = loc2degrees (Location (sembmaxX, sembmaxY), origin)
        azi   = toAzimuth (float(Origin['lat']), float(Origin['lon']),float(sembmaxX), float(sembmaxY))

        sembmaxvaluev[a] = sembmax
        sembmaxlatv[a]   = sembmaxX
        sembmaxlonv[a]   = sembmaxY

        fobjsembmax.write ('%d %.2f %.2f %.20f %.20f %d %03f %f %03f\n' % (a*step,sembmaxX,sembmaxY,sembmax,uncert,usedarrays,delta,float(azi),delta*119.19))
        fobj.close()


    fobjsembmax.close()

    durationpath  = os.path.join (folder, "duration.txt")
    trigger.writeSembMaxValue (sembmaxvaluev,sembmaxlatv,sembmaxlonv,ntimes,Config,Folder)
    trigger.semblancestalta (sembmaxvaluev,sembmaxlatv,sembmaxlonv)

def toMatrix (npVector, nColumns) :

    t   = npVector.tolist()[0]
    n   = nColumns
    mat = []

    for i in range (len(t) / n) :
       pos1  = i * n
       pos2  = pos1 + n
       slice = t [pos1:pos2]
       assert len(slice) == nColumns
       mat.append (slice)

    return mat




def  doCalc (flag,Config,WaveformDict,FilterMetaData,Gmint,Gmaxt,TTTGridMap,Folder,Origin, ntimes, switch, ev,arrayfolder, syn_in):
    '''
    method for calculating semblance of one station array
    '''
    Logfile.add ('PROCESS %d %s' % (flag,' Enters Semblance Calculation') )
    Logfile.add ('MINT  : %f  MAXT: %f Traveltime' % (Gmint,Gmaxt))

    cfg = ConfigObj (dict=Config)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    new_frequence   = cfg.newFrequency()          #('new_frequence')
    forerun= cfg.Int('forerun')
    duration= cfg.Int('duration')
    gridspacing = cfg.Float('gridspacing')

    nostat = len (WaveformDict)
    traveltimes = {}
    recordstarttime = ''
    minSampleCount  = 999999999

    if cfg.UInt ('forerun')>0:
        ntimes = int ((cfg.UInt ('forerun') + cfg.UInt ('duration') ) / cfg.UInt ('step') )
    else:
        ntimes = int ((cfg.UInt ('duration') ) / cfg.UInt ('step') )
    nsamp  = int (winlen * new_frequence)
    nstep  = int (step   * new_frequence)
    from pyrocko import obspy_compat
    from pyrocko import orthodrome, model
    obspy_compat.plant()

    ############################################################################
    calcStreamMap = WaveformDict

    stations = []
    py_trs = []
    for trace in calcStreamMap.keys():
        py_tr = obspy_compat.to_pyrocko_trace(calcStreamMap[trace])
        py_trs.append(py_tr)
        for il in FilterMetaData:
            if str(il) == str(trace):
                        szo = model.Station(lat=il.lat, lon=il.lon,
                                            station=il.sta, network=il.net,
                                            channels=py_tr.channel,
                                            elevation=il.ele, location=il.loc)
                        stations.append(szo) #right number of stations?


#==================================synthetic BeamForming=======================================

    if cfg.Bool('shift_by_phase_pws') == True:
        calcStreamMapshifted= calcStreamMap.copy()
        from obspy.core import stream
        stream = stream.Stream()
        for trace in calcStreamMapshifted.keys():
            stream.append(calcStreamMapshifted[trace])
        pws_stack = PWS_stack([stream], weight=2, normalize=True)
        for tr in pws_stack:
            for trace in calcStreamMapshifted.keys():
                    calcStreamMapshifted[trace]=tr
        calcStreamMap = calcStreamMapshifted


    if cfg.Bool('shift_by_phase_onset') == True:
        pjoin = os.path.join
        timeev = util.str_to_time(ev.time)
        trs_orgs= []
        calcStreamMapshifted= calcStreamMap.copy()
        for trace in calcStreamMapshifted.keys():
                tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapshifted[trace])
                trs_orgs.append(tr_org)

        timing = CakeTiming(
           phase_selection='first(p|P|PP|P(cmb)P(icb)P(icb)p(cmb)p)-20',
           fallback_time=100.)
        traces = trs_orgs

        event = model.Event(lat=float(ev.lat), lon=float(ev.lon), depth=ev.depth*1000., time=timeev)
        directory = arrayfolder
        bf = BeamForming(stations, traces, normalize=True)
        shifted_traces = bf.process(event=event,
                  timing=timing,
                  fn_dump_center=pjoin(directory, 'array_center.pf'),
                  fn_beam=pjoin(directory, 'beam.mseed'))
        i = 0
        store_id = syn_in.store()
        engine = LocalEngine(store_superdirs=[syn_in.store_superdirs()])
        for trace in calcStreamMapshifted.keys():
            recordstarttime = calcStreamMapshifted[trace].stats.starttime.timestamp
            recordendtime = calcStreamMapshifted[trace].stats.endtime.timestamp
            mod = shifted_traces[i]
            extracted = mod.chop(recordstarttime, recordendtime, inplace=False)
            shifted_obs_tr = obspy_compat.to_obspy_trace(extracted)
            calcStreamMapshifted[trace]=shifted_obs_tr
            i = i+1

        calcStreamMap = calcStreamMapshifted


    weight = 0.
    if cfg.Bool('weight_by_noise') == True:
        from noise_analyser import analyse
        pjoin = os.path.join
        timeev = util.str_to_time(ev.time)
        trs_orgs= []
        calcStreamMapshifted= calcStreamMap.copy()
        for trace in calcStreamMapshifted.keys():
                tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapshifted[trace])
                trs_orgs.append(tr_org)

        timing = CakeTiming(
           phase_selection='first(p|P|PP|P(cmb)P(icb)P(icb)p(cmb)p)-20',
           fallback_time=100.)
        traces = trs_orgs
        event = model.Event(lat=float(ev.lat), lon=float(ev.lon), depth=ev.depth*1000., time=timeev)
        directory = arrayfolder
        bf = BeamForming(stations, traces, normalize=True)
        shifted_traces = bf.process(event=event,
                  timing=timing,
                  fn_dump_center=pjoin(directory, 'array_center.pf'),
                  fn_beam=pjoin(directory, 'beam.mseed'))
        i = 0
        store_id = syn_in.store()
        engine = LocalEngine(store_superdirs=[syn_in.store_superdirs()])
        weight = analyse(shifted_traces, engine, event, stations,
         100., store_id, nwindows=1,
         check_events=True, phase_def='P')

    for trace in calcStreamMap.keys():
        recordstarttime = calcStreamMap[trace].stats.starttime
        d = calcStreamMap[trace].stats.starttime
        d = d.timestamp

        if calcStreamMap[trace].stats.npts < minSampleCount:
            minSampleCount = calcStreamMap[trace].stats.npts

    ############################################################################
    traces = num.ndarray (shape=(len(calcStreamMap), minSampleCount), dtype=float)
    traveltime = num.ndarray (shape=(len(calcStreamMap), dimX*dimY), dtype=float)
    latv   = num.ndarray (dimX*dimY, dtype=float)
    lonv   = num.ndarray (dimX*dimY, dtype=float)
    ############################################################################


    c=0
    streamCounter = 0

    for key in calcStreamMap.keys():
        streamID = key
        c2   = 0

        for o in calcStreamMap[key]:
            if c2 < minSampleCount:
                traces[c][c2] = o

                c2 += 1


        for key in TTTGridMap.keys():

            if streamID == key:
                traveltimes[streamCounter] = TTTGridMap[key]
            else:
                "NEIN", streamID, key


        if not streamCounter in traveltimes :
           continue                              #hs : thread crashed before

        g = traveltimes[streamCounter]
        dimZ  = g.dimZ
        mint  = g.mint
        maxt  = g.maxt
        Latul = g.Latul
        Lonul = g.Lonul
        Lator = g.Lator
        Lonor = g.Lonor

        gridElem = g.GridArray

        for x in range(dimX):
            for y in range(dimY):
                elem = gridElem[x, y]

                traveltime [c][x * dimY + y] = elem.tt
                latv [x * dimY + y] = elem.lat
                lonv [x * dimY + y] = elem.lon
        #endfor

        c += 1
        streamCounter += 1

    #endfor


# ==================================semblance calculation=======

    t1 = time.time()
    traces = traces.reshape(1, nostat*minSampleCount)

    traveltimes = traveltime.reshape(1, nostat*dimX*dimY)
    TTTGrid = True
    manual_shift = False

    if manual_shift:

        pjoin = os.path.join
        timeev = util.str_to_time(ev.time)
        trs_orgs = []
        calcStreamMapshifted = calcStreamMap.copy()
        for trace in calcStreamMapshifted.keys():
                tr_org = obspy_compat.to_pyrocko_trace(
                    calcStreamMapshifted[trace])
                trs_orgs.append(tr_org)

        timing = CakeTiming(
           phase_selection='first(p|P|PP|P(cmb)P(icb)P(icb)p(cmb)p)-20',
           fallback_time=100.)
        traces = trs_orgs
        backSemb = num.ndarray(shape=(ntimes, dimX*dimY), dtype=float)
        bf = BeamForming(stations, traces, normalize=True)

        for i in range(ntimes):
            sembmax = 0
            sembmaxX = 0
            sembmaxY = 0
            for j in range(dimX * dimY):
                event = model.Event(lat=float(latv[j]), lon=float(lonv[j]),
                                    depth=ev.depth*1000., time=timeev)
                directory = arrayfolder
                shifted_traces, stack = bf.process(event=event,
                                                   timing=timing,
                                                   fn_dump_center=pjoin(
                                                                directory,
                                                         'array_center.pf'),
                                                   fn_beam=pjoin(directory,
                                                                 'beam.mseed'))
                tmin = stack.tmin+(i*nstep)+20
                tmax = stack.tmin+(i*nstep)+60
                stack.chop(tmin, tmax)
                backSemb[i][j] = abs(sum(stack.ydata))

        k = backSemb
        TTTGrid = False

    if TTTGrid:
        start_time = time.time()
        if cfg.UInt('forerun') > 0:
            ntimes = int((cfg.UInt('forerun') + cfg.UInt('duration'))/step)
        else:
            ntimes = int((cfg.UInt('duration')) / step)
        nsamp = int(winlen)
        nstep = int(step)
        Gmint = cfg.Int('forerun')

        k = semblance(maxp, nostat, nsamp, ntimes, nstep, dimX, dimY, Gmint,
                      new_frequence, minSampleCount, latv, lonv, traveltimes,
                      traces, calcStreamMap, timeev, Config, Origin)
        print("--- %s seconds ---" % (time.time() - start_time))

    t2 = time.time()

    Logfile.add('%s took %0.3f s' % ('CALC:',(t2-t1)))

    partSemb = k
    partSemb = partSemb.reshape(ntimes, migpoints)

    return partSemb


def  doCalc_syn (flag,Config,WaveformDict,FilterMetaData,Gmint,Gmaxt,TTTGridMap,
                Folder,Origin, ntimes, switch, ev,arrayfolder, syn_in, parameter):
    '''
    method for calculating semblance of one station array
    '''
    Logfile.add ('PROCESS %d %s' % (flag,' Enters Semblance Calculation') )
    Logfile.add ('MINT  : %f  MAXT: %f Traveltime' % (Gmint,Gmaxt))

    cfg = ConfigObj (dict=Config)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    new_frequence   = cfg.newFrequency()          #('new_frequence')
    forerun= cfg.Int('forerun')
    duration= cfg.Int('duration')
    gridspacing = cfg.Float('gridspacing')

    nostat = len (WaveformDict)
    traveltimes = {}
    recordstarttime = ''
    minSampleCount  = 999999999

    if cfg.UInt ('forerun')>0:
        ntimes = int ((cfg.UInt ('forerun') + cfg.UInt ('duration') ) / cfg.UInt ('step') )
    else:
        ntimes = int ((cfg.UInt ('duration') ) / cfg.UInt ('step') )
    nsamp  = int (winlen * new_frequence)
    nstep  = int (step   * new_frequence)
    from pyrocko import obspy_compat
    from pyrocko import orthodrome, model
    obspy_compat.plant()

    ############################################################################
    calcStreamMap = WaveformDict

    stations = []
    py_trs = []
    for trace in calcStreamMap.keys():
        py_tr = obspy_compat.to_pyrocko_trace(calcStreamMap[trace])
        py_trs.append(py_tr)
        for il in FilterMetaData:
            if str(il) == str(trace):
                        szo = model.Station(lat=il.lat, lon=il.lon,
                                            station=il.sta, network=il.net,
                                            channels=py_tr.channel,
                                            elevation=il.ele, location=il.loc)
                        stations.append(szo) #right number of stations?

    store_id = syn_in.store()
    engine = LocalEngine(store_superdirs=[syn_in.store_superdirs()])

    targets = []
    for st in stations:
        target = Target(
                lat=st.lat,
                lon=st.lon,
                store_id=store_id,
                codes=(st.network, st.station, st.location, 'BHZ'),
                tmin=-1900,
                tmax=3900,
                interpolation='multilinear',
                quantity=cfg.quantity())
        targets.append(target)

    if syn_in.nsources() == 1:
        if syn_in.use_specific_stf() is True:
            stf = syn_in.stf()
            exec(stf)
        else:
            stf = STF()
        if syn_in.source() == 'RectangularSource':
                source = RectangularSource(
                    lat=float(syn_in.lat_0()),
                    lon=float(syn_in.lon_0()),
                    depth=syn_in.depth_syn_0()*1000.,
                    strike=syn_in.strike_0(),
                    dip=syn_in.dip_0(),
                    rake=syn_in.rake_0(),
                    width=syn_in.width_0()*1000.,
                    length=syn_in.length_0()*1000.,
                    nucleation_x=syn_in.nucleation_x_0(),
                    slip=syn_in.slip_0(),
                    nucleation_y=syn_in.nucleation_y_0(),
                    stf=stf,
                    time=util.str_to_time(syn_in.time_0()))
        if syn_in.source() == 'DCSource':
                source = DCSource(
                    lat=float(syn_in.lat_0()),
                    lon=float(syn_in.lon_0()),
                    depth=syn_in.depth_syn_0()*1000.,
                    strike=syn_in.strike_0(),
                    dip=syn_in.dip_0(),
                    rake=syn_in.rake_0(),
                    stf=stf,
                    time=util.str_to_time(syn_in.time_0()),
                    magnitude=syn_in.magnitude_0())

    else:
        sources = []
        for i in range(syn_in.nsources()):
            if syn_in.use_specific_stf() is True:
                stf = syn_in.stf()
                exec(stf)

            else:
                stf = STF()
            if syn_in.source() == 'RectangularSource':
                    sources.append(RectangularSource(
                        lat=float(syn_in.lat_1(i)),
                        lon=float(syn_in.lon_1(i)),
                        depth=syn_in.depth_syn_1(i)*1000.,
                        strike=syn_in.strike_1(i),
                        dip=syn_in.dip_1(i),
                        rake=syn_in.rake_1(i),
                        width=syn_in.width_1(i)*1000.,
                        length=syn_in.length_1(i)*1000.,
                        nucleation_x=syn_in.nucleation_x_1(i),
                        slip=syn_in.slip_1(i),
                        nucleation_y=syn_in.nucleation_y_1(i),
                        stf=stf,
                        time=util.str_to_time(syn_in.time_1(i))))

            if syn_in.source() == 'DCSource':
                    sources.append(DCSource(
                        lat=float(syn_in.lat_1(i)),
                        lon=float(syn_in.lon_1(i)),
                        depth=syn_in.depth_1(i)*1000.,
                        strike=syn_in.strike_1(i),
                        dip=syn_in.dip_1(i),
                        rake=syn_in.rake_1(i),
                        stf=stf,
                        time=util.str_to_time(syn_in.time_1(i)),
                        magnitude=syn_in.magnitude_1(i)))
        source = CombiSource(subsources=sources)
    response = engine.process(source, targets)

    synthetic_traces = response.pyrocko_traces()
    if cfg.Bool('synthetic_test_add_noise') is True:
        from noise_addition import add_noise
        trs_orgs = []
        calcStreamMapsyn = calcStreamMap.copy()
        #from pyrocko import trace
        for tracex in calcStreamMapsyn.keys():
                for trl in synthetic_traces:
                    if str(trl.name()[4:12]) == str(tracex[4:]):
                        tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapsyn[tracex])
                        tr_org.downsample_to(2.0)
                        trs_orgs.append(tr_org)
        store_id = syn_in.store()
        engine = LocalEngine(store_superdirs=[syn_in.store_superdirs()])
        synthetic_traces = add_noise(trs_orgs, engine, source.pyrocko_event(),
                                     stations,
                                     store_id, phase_def='P')
    trs_org = []
    trs_orgs = []
    fobj = os.path.join(arrayfolder, 'shift.dat')
    xy = num.loadtxt(fobj, usecols=1, delimiter=',')
    calcStreamMapsyn = calcStreamMap.copy()
    #from pyrocko import trace
    for tracex in calcStreamMapsyn.keys():
            for trl in synthetic_traces:
                if str(trl.name()[4:12])== str(tracex[4:]):
                    mod = trl

                    recordstarttime = calcStreamMapsyn[tracex].stats.starttime.timestamp
                    recordendtime = calcStreamMapsyn[tracex].stats.endtime.timestamp
                    tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapsyn[tracex])
                    trs_orgs.append(tr_org)

                    tr_org_add = mod.chop(recordstarttime, recordendtime, inplace=False)
                    synthetic_obs_tr = obspy_compat.to_obspy_trace(tr_org_add)
                    calcStreamMapsyn[tracex] = synthetic_obs_tr
                    trs_org.append(tr_org_add)
    calcStreamMap = calcStreamMapsyn

    if cfg.Bool('shift_by_phase_pws') == True:
        calcStreamMapshifted= calcStreamMap.copy()
        from obspy.core import stream
        stream = stream.Stream()
        for trace in calcStreamMapshifted.keys():
            stream.append(calcStreamMapshifted[trace])
        pws_stack = PWS_stack([stream], weight=2, normalize=True)
        for tr in pws_stack:
            for trace in calcStreamMapshifted.keys():
                    calcStreamMapshifted[trace]=tr
        calcStreamMap = calcStreamMapshifted


    if cfg.Bool('shift_by_phase_onset') == True:
        pjoin = os.path.join
        timeev = util.str_to_time(ev.time)
        trs_orgs= []
        calcStreamMapshifted= calcStreamMap.copy()
        for trace in calcStreamMapshifted.keys():
                tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapshifted[trace])
                trs_orgs.append(tr_org)

        timing = CakeTiming(
           phase_selection='first(p|P|PP|P(cmb)P(icb)P(icb)p(cmb)p)-20',
           fallback_time=100.)
        traces = trs_orgs

        event = model.Event(lat=float(ev.lat), lon=float(ev.lon), depth=ev.depth*1000., time=timeev)
        directory = arrayfolder
        bf = BeamForming(stations, traces, normalize=True)
        shifted_traces = bf.process(event=event,
                  timing=timing,
                  fn_dump_center=pjoin(directory, 'array_center.pf'),
                  fn_beam=pjoin(directory, 'beam.mseed'))
        i = 0
        store_id = syn_in.store()
        engine = LocalEngine(store_superdirs=[syn_in.store_superdirs()])
        for trace in calcStreamMapshifted.keys():
            recordstarttime = calcStreamMapshifted[trace].stats.starttime.timestamp
            recordendtime = calcStreamMapshifted[trace].stats.endtime.timestamp
            mod = shifted_traces[i]
            extracted = mod.chop(recordstarttime, recordendtime, inplace=False)
            shifted_obs_tr = obspy_compat.to_obspy_trace(extracted)
            calcStreamMapshifted[trace]=shifted_obs_tr
            i = i+1

        calcStreamMap = calcStreamMapshifted


    weight = 0.
    if cfg.Bool('weight_by_noise') == True:
        from noise_analyser import analyse
        pjoin = os.path.join
        timeev = util.str_to_time(ev.time)
        trs_orgs= []
        calcStreamMapshifted= calcStreamMap.copy()
        for trace in calcStreamMapshifted.keys():
                tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapshifted[trace])
                trs_orgs.append(tr_org)

        timing = CakeTiming(
           phase_selection='first(p|P|PP|P(cmb)P(icb)P(icb)p(cmb)p)-20',
           fallback_time=100.)
        traces = trs_orgs
        event = model.Event(lat=float(ev.lat), lon=float(ev.lon), depth=ev.depth*1000., time=timeev)
        directory = arrayfolder
        bf = BeamForming(stations, traces, normalize=True)
        shifted_traces = bf.process(event=event,
                  timing=timing,
                  fn_dump_center=pjoin(directory, 'array_center.pf'),
                  fn_beam=pjoin(directory, 'beam.mseed'))
        i = 0
        store_id = syn_in.store()
        engine = LocalEngine(store_superdirs=[syn_in.store_superdirs()])
        weight = analyse(shifted_traces, engine, event, stations,
         100., store_id, nwindows=1,
         check_events=True, phase_def='P')

    for trace in calcStreamMap.keys():
        recordstarttime = calcStreamMap[trace].stats.starttime
        d = calcStreamMap[trace].stats.starttime
        d = d.timestamp

        if calcStreamMap[trace].stats.npts < minSampleCount:
            minSampleCount = calcStreamMap[trace].stats.npts

    ############################################################################
    traces = num.ndarray (shape=(len(calcStreamMap), minSampleCount), dtype=float)
    traveltime = num.ndarray (shape=(len(calcStreamMap), dimX*dimY), dtype=float)
    latv   = num.ndarray (dimX*dimY, dtype=float)
    lonv   = num.ndarray (dimX*dimY, dtype=float)
    ############################################################################


    c=0
    streamCounter = 0

    for key in calcStreamMap.keys():
        streamID = key
        c2   = 0

        for o in calcStreamMap[key]:
            if c2 < minSampleCount:
                traces[c][c2] = o

                c2 += 1


        for key in TTTGridMap.keys():

            if streamID == key:
                traveltimes[streamCounter] = TTTGridMap[key]
            else:
                "NEIN", streamID, key


        if not streamCounter in traveltimes :
           continue                              #hs : thread crashed before

        g = traveltimes[streamCounter]
        dimZ  = g.dimZ
        mint  = g.mint
        maxt  = g.maxt
        Latul = g.Latul
        Lonul = g.Lonul
        Lator = g.Lator
        Lonor = g.Lonor

        gridElem = g.GridArray

        for x in range(dimX):
            for y in range(dimY):
                elem = gridElem[x, y]

                traveltime [c][x * dimY + y] = elem.tt
                latv [x * dimY + y] = elem.lat
                lonv [x * dimY + y] = elem.lon
        #endfor

        c += 1
        streamCounter += 1

    #endfor


    ############################## CALCULATE PARAMETER FOR SEMBLANCE CALCULATION ##################
    nsamp = winlen * new_frequence

    nstep = int (step*new_frequence)
    migpoints = dimX * dimY

    dimZ = 0
    new_frequence = cfg.newFrequency ()              # ['new_frequence']
    maxp = int (Config['ncore'])


    Logfile.add ('PROCESS %d  NTIMES: %d' % (flag,ntimes))

    if False :
       print ('nostat ',nostat,type(nostat))
       print ('nsamp ',nsamp,type(nsamp))
       print ('ntimes ',ntimes,type(ntimes))
       print ('nstep ',nstep,type(nstep))
       print ('dimX ',dimX,type(dimX))
       print ('dimY ',dimY,type(dimY))
       print ('mint ',Gmint,type(mint))
       print ('new_freq ',new_frequence,type(new_frequence))
       print ('minSampleCount ',minSampleCount,type(minSampleCount))
       print ('latv ',latv,type(latv))
       print ('traces',traces,type(traces))
       print ('traveltime',traveltime,type(traveltime))


#==================================semblance calculation========================================

    t1 = time.time()
    traces = traces.reshape   (1,nostat*minSampleCount)
    traveltime = traveltime.reshape (1,nostat*dimX*dimY)
    USE_C_CODE = True
    try:
        if USE_C_CODE :
            import Cm
            import CTrig
            start_time = time.time()
            k  = Cm.otest (maxp,nostat,nsamp,ntimes,nstep,dimX,dimY,Gmint,new_frequence,
                          minSampleCount,latv,lonv,traveltime,traces)
            print("--- %s seconds ---" % (time.time() - start_time))
        else :
            start_time = time.time()
            k = otest (maxp,nostat,nsamp,ntimes,nstep,dimX,dimY,Gmint,new_frequence,
                      minSampleCount,latv,lonv,traveltime,traces)                       #hs
            print("--- %s seconds ---" % (time.time() - start_time))
    except:
        print("loaded tttgrid has probably wrong dimensions or stations, delete\
                ttgrid or exchange")

    t2 = time.time()


    partSemb = k

    partSemb_syn  = partSemb.reshape (ntimes,migpoints)


    return partSemb_syn

def optimization(*params, **args):
    counter = params[1]
    Config = params[2]
    Wdf = params[3]
    FilterMeta = params[4]
    mint = params[5]
    maxt = params[6]
    TTTGridMap = params[7]
    Folder = params[8]
    Origin = params[9]
    ntimes = params[10]
    switch = params[11]
    ev = params[12]
    arrayfolder = params[13]
    syn_in = params[14]
    data = params[15]
    evpath = params[16]
    XDict = params[17]
    RefDict = params[18]
    workdepth = params[19]
    filterindex = params[20]
    Wdfs = params[21]

    networks = Config['networks'].split(',')
    params = num.asarray(params)
    parameter = num.ndarray.tolist(params)
    ASL_syn = []


    C  = config.Config (evpath)
    Config = C.parseConfig ('config')
    cfg = ConfigObj (dict=Config)
    if cfg.pyrocko_download() == True:
        Meta = C.readpyrockostations()#

    elif cfg.colesseo_input() == True:
        scenario = guts.load(filename=cfg.colosseo_scenario_yml())
        scenario_path = cfg.colosseo_scenario_yml()[:-12]
        Meta = C.readcolosseostations(scenario_path)
    else:
        Meta = C.readMetaInfoFile()
    l = 0
    for i in networks:

        arrayname = i
        arrayfolder = os.path.join (Folder['semb'],arrayname)

        network = Config[i].split('|')

        FilterMeta = ttt.filterStations (Meta,Config,Origin,network)

        if len(FilterMeta)  < 3: continue

        W = XDict[i]
        refshift = RefDict[i]

        FilterMeta = cmpFilterMetavsXCORR (W, FilterMeta)

        Logfile.add ('BOUNDING BOX DIMX: %s  DIMY: %s  GRIDSPACING: %s \n'
                 % (Config['dimx'],Config['dimy'],Config['gridspacing']))

        f = open('../tttgrid/tttgrid_%s_%s_%s.pkl' % (ev.time, arrayname, workdepth), 'rb')
        TTTGridMap,mint,maxt = pickle.load(f)
        f.close()


        switch = filterindex

        tw  = times.calculateTimeWindows (mint,maxt,Config,ev, switch)
        Wdf = Wdfs[l]
        semb_syn = doCalc_syn (counter,Config,Wdf,FilterMeta,mint,maxt,TTTGridMap,
                                     Folder,Origin,ntimes,switch, ev,arrayfolder, syn_in,
                                      parameter[0])
        ASL_syn.append(semb_syn)
        counter += 1
        l += 1

    sembmax_syn = sembCalc.collectSemb(ASL_syn,Config,Origin,Folder,ntimes,len(networks),switch)

    misfit_list = []  # init a list for a all the singular misfits
    norm_list = []  # init a list for a all the singular normalizations
    taper = trace.CosFader(xfade=2.0)  # Cosine taper with fade in and out of 2s.
    bw_filter = trace.ButterworthResponse(corner=0.000055,  # in Hz
                                      order=4,
                                      type='high')  # "low"pass or "high"pass
    setup = trace.MisfitSetup(description='Misfit Setup',
                              norm=2,  # L1 or L2 norm
                              taper=taper,
                              filter=bw_filter,
                              domain='time_domain')
    nsamples = len(data)
    tmin = util.str_to_time('2010-02-20 15:15:30.100')
    tr = trace.Trace(station='TEST', channel='Z',
                     deltat=0.5, tmin=tmin, ydata=data)
    syn = trace.Trace(station='TEST', channel='Z',
                     deltat=0.5, tmin=tmin, ydata=sembmax_syn)
    misfit, norm = tr.misfit(candidate=syn, setup=setup) # calculate the misfit of a single observed trace with its synthetics
    # with the setup from above
    misfit_list.append(misfit), norm_list.append(norm)  # append the misfit into a list
    global_misfit_normed = num.sqrt(num.nansum((num.asarray(misfit_list))**2) / # sum all the misfits and normalize to get a single minimizable value
                                    num.nansum((num.asarray(norm_list))**2))
    return global_misfit_normed


def solve(counter,Config,Wdf,FilterMeta,mint,maxt,TTTGridMap,
                             Folder,Origin,ntimes,switch, ev,arrayfolder,
                             syn_in, ASL_d, sembmax_d, evpath, XDict,
                             RefDict, workdepth, filterindex, Wdfs):
    import scipy
    t = time.time()  # start timing
    # bounds given as (min,max)

    bounds = ((syn_in.mag_0_low(), syn_in.mag_0_high()),  # magnitude
            (syn_in.strike_0_low(), syn_in.strike_0_high()),  # strike [deg.]
            (syn_in.dip_0_low(), syn_in.dip_0_high()),  # dip [deg.]
            (syn_in.rake_0_low(), syn_in.rake_0_high()),  # rake [deg.]
            (syn_in.depth_0_low()*km, syn_in.depth_0_high()*km),  # depth [km]
            (syn_in.north_shift_0_low()*km, syn_in.north_shift_0_high()*km),  # north shift from GCMT [km]
            (syn_in.east_shift_0_low()*km, syn_in.east_shift_0_high()*km),  # east shift from GCMT [km]
            (syn_in.time_0_low(), syn_in.time_0_high()))  # timeshift from GCMT [s]
    # optimize.differential_evolution of scipy is used for the optim.
    # Differential Evolution is stochastic in nature (does not use gradient methods)
    #to find the minimium, and can search large areas of candidate space, but often requires
    #larger numbers of function evaluations than conventional gradient based techniques.
    # The scipy solver can easily be exchanged.

    data = sembmax_d
    result = scipy.optimize.differential_evolution(optimization, bounds=bounds, maxiter=3, popsize=3, args=(counter,Config,Wdf,FilterMeta,mint,maxt,TTTGridMap,
                                 Folder,Origin,ntimes,switch, ev,arrayfolder, syn_in, data, evpath, XDict, RefDict, workdepth, filterindex, Wdfs))
    elapsed = time.time() - t  # get the processing time
    print(result.x)
