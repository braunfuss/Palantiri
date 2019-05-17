import os
import sys
import math
from math import radians, cos, sin, atan2
import logging
sys.path.append('../tools/')
sys.path.append('../Common/')
from obspy.core.utcdatetime import UTCDateTime
import Basic
import Globals
import Logfile
import DataTypes
from DataTypes import Location
from ObspyFkt import loc2degrees
from pyrocko import orthodrome, obspy_compat, model
from ConfigFile import ConfigObj, OriginCfg, SynthCfg, FilterCfg
import time
import numpy as num
from collections import defaultdict
from pyrocko.gf import ws, LocalEngine, Target, DCSource, RectangularSource, MTSource
from pyrocko import util, pile, model, catalog, gf, cake
from pyrocko.guts import Object, String, Float, List
km = 1000.
r2d = 180./math.pi
d2r = 1./r2d
earthradius = 6371.*km
import trigger
from semp import semblance
from beam_stack import BeamForming
from pyrocko.gf import STF
from stacking import PWS_stack
import copy
if sys.version_info.major >= 3:
    xrange = range

logger = logging.getLogger('ARRAY-MP')


def make_bayesian_weights(narrays, nbootstrap=100,
                          type='bayesian', rstate=None):

    ws = num.zeros((int(nbootstrap), int(narrays)))
    if rstate is None:
        rstate = num.random.RandomState()

    for ibootstrap in range(nbootstrap):
        if type == 'classic':
            ii = rstate.randint(0, narrays, size=narrays)
            ws[ibootstrap, :] = num.histogram(
                ii, narrays, (-0.5, narrays - 0.5))[0]
        elif type == 'bayesian':
            f = rstate.uniform(0., 1., size=narrays+1)
            f[0] = 0.
            f[-1] = 1.
            f = num.sort(f)
            g = f[1:] - f[:-1]
            ws[ibootstrap, :] = g
        else:
            assert False
    return ws


def spectrum(data, deltat, pad_to_pow2=False, tfade=None):
        '''
        Get FFT spectrum of trace.
        :param pad_to_pow2: whether to zero-pad the data to next larger
            power-of-two length
        :param tfade: ``None`` or a time length in seconds, to apply cosine
            shaped tapers to both
        :returns: a tuple with (frequencies, values)
        '''

        ndata = data.size

        if pad_to_pow2:
            ntrans = nextpow2(ndata)
        else:
            ntrans = ndata

        if tfade is None:
            ydata = data
        else:
            ydata = self.ydata * costaper(
                0., tfade, deltat*(ndata-1)-tfade, deltat*ndata,
                ndata, deltat)

        fydata = num.fft.fft(ydata)
        df = 1./(ntrans*deltat)
        fxdata = num.arange(len(fydata))*df
        return fxdata, fydata


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
        :param z_dist: tuple with(depth, distance)
        '''
        z, dist = z_dist
        if(dist, z) in self.arrivals.keys():
            return self.return_time(self.arrivals[(dist, z)])

        phases = [cake.PhaseDef(pid) for pid in self.phases]
        arrivals = mod.arrivals(
            distances=[dist*cake.m2d], phases=phases, zstart=z)
        if arrivals == []:
            logger.warn(
                'no phase at d=%s, z=%s.(return fallback time)' %(dist, z))
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


class SembMax(object):
    '''
    class to store sembmax object for each grid point
    '''

    def __init__(self, lat, lon, semb):

        self.lat = lat
        self.lon = lon
        self.semb = semb


class FileSembMax(object):
    '''
    class to strore sembmax object for the sembmaxvalue file
    '''

    def __init__(self,istep,sembmaxX,sembmaxY,sembmax,usedarrays,delta,azi,deltakm):

        self.istep = istep
        self.sembmaxX = sembmaxX
        self.sembmaxY = sembmaxY
        self.sembmax = sembmax
        self.usedarrays = usedarrays
        self.delta = delta
        self.azi= azi
        self.deltakm= deltakm

    def get(self):
        return('%d %.2f %.2f %f %d %03f %f %03f\n' %(self.istep,self.sembmaxX,self.sembmaxY,
                                                       self.sembmax,self.usedarrays,self.delta,
                                                       self.azi,self.delta*119.19))


def toAzimuth(latevent, lonevent, latsource, lonsource):
        '''
        method to calculate azimuth between two points
        '''

        lat1 = radians(latsource)
        lon1 = radians(lonsource)
        lat2 = radians(latevent)
        lon2 = radians(lonevent)

        x = sin(lon1-lon2 ) * cos(lat2)
        y = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1-lon2)
        angle = -atan2(x, y)

        if angle < 0.0:
         angle  += math.pi * 2.0

        angle = math.degrees(angle)
        angle = '%02f' % angle

        return angle


def writeSembMatricesSingleArray(SembList, Config, Origin, arrayfolder, ntimes,
                                 switch, phase, bootstrap=None):
    '''
    method to write semblance matrizes from one processes to file for each
    timestep.
    '''
    logger.info('start write semblance matrices')

    cfg = ConfigObj(dict=Config)
    origin = OriginCfg(Origin)

    dimX = cfg.dimX()
    dimY = cfg.dimY()
    if switch == 0:
        winlen = cfg.Float('winlen')
        step = cfg.Float('step')
    if switch == 1:
        winlen = cfg.Float('winlen_f2')
        step = cfg.Float('step_f2')

    latv = []
    lonv = []

    gridspacing = cfg.Float('gridspacing')
    migpoints = dimX * dimY

    o_lat = origin.lat()
    o_lon = origin.lon()
    oLatul = 0
    oLonul = 0

    z = 0

    for i in xrange(dimX):
        oLatul = o_lat - ((dimX-1)/2) * gridspacing + i*gridspacing

        if z == 0 and i == 0:
            Latul = oLatul
        o = 0

        for j in xrange(dimY):
            oLonul = o_lon - ((dimY-1)/2) * gridspacing + j * gridspacing

            if o == 0 and j == 0:
                Lonul = oLonul

            latv.append(oLatul)
            lonv.append(oLonul)

    rc = UTCDateTime(Origin['time'])
    rcs = '%s-%s-%s_%02d:%02d:%02d' % (rc.day, rc.month, rc.year, rc.hour,
                                       rc.minute, rc.second)
    d = rc.timestamp

    for a, i in enumerate(SembList):
        if bootstrap is None:
            fobj = open(os.path.join(arrayfolder, '%s-%s_%03d_%s.ASC'
                                     % (switch, Origin['depth'], a,
                                        phase)), 'w')
        else:
            fobj = open(os.path.join(arrayfolder, '%s-%s-boot%s_%03d_%s.ASC'
                                     % (switch, Origin['depth'], bootstrap,
                                        a, phase)), 'w')

        fobj.write('# %s , %s\n' % (d, rcs))
        fobj.write('# step %ds| ntimes %d| winlen: %ds\n' % (step, ntimes,
                                                             winlen))
        fobj.write('# \n')
        fobj.write('# southwestlat: %.2f dlat: %f nlat: %f \n' % (Latul,
                                                                  gridspacing,
                                                                  dimX))
        fobj.write('# southwestlon: %.2f dlon: %f nlon: %f \n' % (Lonul,
                                                                  gridspacing,
                                                                  dimY))
        fobj.write('# ddepth: 0 ndepth: 1 \n')

        for j in range(migpoints):
            x = latv[j]
            y = lonv[j]
            z = origin.depth()
            semb = i[j]

            fobj.write('%.2f %.2f %.2f %.20f\n' % (x, y, z, semb))

        fobj.close()



def collectSemb(SembList, Config, Origin, Folder, ntimes, arrays, switch,
                array_centers, phase, cboot=None, temp_comb=None):
    '''
    method to collect semblance matrizes from all processes and write them to file for each timestep
    '''
    Logfile.add('start collect in collectSemb')
    cfg = ConfigObj(dict=Config)
    origin = ConfigObj(dict=Origin)

    dimX = cfg.dimX()
    dimY = cfg.dimY()
    if switch == 0:
        winlen = cfg.winlen()
        step = cfg.step()
    if switch == 1:
        winlen = cfg.winlen_f2()
        step = cfg.step_f2()

    latv = []
    lonv = []

    gridspacing = cfg.Float('gridspacing')
    migpoints = dimX * dimY
    o_lat = origin.lat()
    o_lon = origin.lon()
    oLatul = 0
    oLonul = 0

    z=0

    for i in xrange(dimX):
         oLatul = o_lat -((dimX-1)/2) * gridspacing + i*gridspacing

         if z == 0 and i == 0 :
             Latul = oLatul
         o=0

         for j in xrange(dimY):
               oLonul = o_lon -((dimY-1)/2) * gridspacing + j*gridspacing

               if o==0 and j==0:
                    Lonul = oLonul

               latv.append(oLatul)
               lonv.append(oLonul)

    origin = DataTypes.dictToLocation(Origin)
    icount = 0

    azis = []
    for a in SembList:
        x = array_centers[icount][0]
        y = array_centers[icount][1]
        delta = orthodrome.distance_accurate50m_numpy(x, y, origin.lat,
                                                      origin.lon)
        azis.append(toAzimuth(float(Origin['lat']), float(Origin['lon']),
                              x, y))
        icount = icount+1
    aziblocks = num.arange(0., 360., 20.)
    aziblockslist = list(aziblocks)
    aziweights = []
    blockweight = 1.
    idxs = []
    tmp_general = 1
    for azi in azis:
            val = min(aziblocks, key=lambda x: abs(x-float(azi)))
            idxs.append(aziblockslist.index(val))
    for idx in idxs:
        aziweights.append(1./idxs.count(idx))

    min_coor = num.zeros([icount, 2])
    if cfg.Bool('bootstrap_array_weights') is True and\
       cfg.Bool('combine_all') is False:
        nboot = cfg.Int('n_bootstrap')
        bs_weights = make_bayesian_weights(len(azis), nbootstrap=nboot)
    else:
        nboot = 1

    if temp_comb is not None:
        nboot = 0
    else:
        nboot = 1

    folder = Folder['semb']
    sembmaxvaluev = num.ndarray(ntimes, dtype=float)
    sembmaxlatv = num.ndarray(ntimes, dtype=float)
    sembmaxlonv = num.ndarray(ntimes, dtype=float)

    rc = UTCDateTime(Origin['time'])
    rcs = '%s-%s-%s_%02d:%02d:%02d' % (rc.day, rc.month, rc.year, rc.hour,
                                       rc.minute, rc.second)
    d = rc.timestamp

    usedarrays = arrays
    for boot in range(nboot):
        c = 0
        tmp = 1

        for a in SembList:
            if cfg.Bool('weight_by_azimuth') is True:
                tmp *= a*aziweights[c]
            elif cfg.Bool('bootstrap_array_weights') is True:
                tmp *= a*bs_weights[boot, c]
            else:
                tmp *= a
            tmp_general *= tmp
            deltas = []
            x = array_centers[c][0]
            y = array_centers[c][1]

            for k in range(0, len(latv)):
                delta = orthodrome.distance_accurate50m_numpy(x, y, latv[k],
                                                              lonv[k])
                deltas.append(orthodrome.distance_accurate50m_numpy(x, y,
                                                                    latv[k],
                                                                    lonv[k]))
                if delta <= num.min(deltas):
                    min_coor[c] = [latv[k], lonv[k]]
            c =+ 1

        array_overlap = num.average(min_coor, axis=0)
        delta_center = orthodrome.distance_accurate50m_numpy(array_overlap[0],
                                                             array_overlap[1],
                                                             origin.lat,
                                                             origin.lon)

        diff_center_lat = origin.lat-array_overlap[0]
        diff_center_lon = origin.lon-array_overlap[1]


        fobjsembmax = open(os.path.join(folder, 'sembmax_%s_boot%s_%s.txt'
                                        % (switch, boot, phase)), 'w')
        norm = num.max(num.max(tmp, axis=1))
        max_p = 0.
        sum_i = 0.

        ishape = tmp[0]
        semb_cum = num.zeros(num.shape(ishape))
        times_cum = num.zeros(num.shape(ishape))
        times_min = num.zeros(num.shape(ishape))
        times_max = num.zeros(num.shape(ishape))
        semb_prior = num.zeros(num.shape(ishape))

    #    correct for array center bias
    #    for j in range(migpoints):
    #                latv[j] = latv[j]#+diff_center_lat
    #                lonv[j] = lonv[j]#+diff_center_lon


        for a, i in enumerate(tmp):
            logger.info('timestep %d' % a)
            fobj = open(os.path.join(folder, '%s-%s_boot%s_%03d_%s.ASC'
                                     % (switch, Origin['depth'], boot, a,
                                        phase)),
                        'w')
            fobj.write('# %s , %s\n' % (d, rcs))
            fobj.write('# step %ds| ntimes %d| winlen: %ds\n'
                       % (step, ntimes, winlen))
            fobj.write('# \n')
            fobj.write('# southwestlat: %.2f dlat: %f nlat: %f \n'
                       % (Latul, gridspacing, dimX))
            fobj.write('# southwestlon: %.2f dlon: %f nlon: %f \n'
                       % (Lonul, gridspacing, dimY))
            fobj.write('# ddepth: 0 ndepth: 1 \n')

            sembmax = 0
            sembmaxX = 0
            sembmaxY = 0
            counter_time = 0
            uncert = num.std(i)
            semb_cum =+ i
            for j in range(migpoints):

                x = latv[j]
                y = lonv[j]

                if cfg.Bool('norm_all') is True:
                    semb = i[j]/norm
                else:
                    semb = i[j]
                fobj.write('%.2f %.2f %.20f\n' % (x, y, semb))
                if semb_prior[j] <= semb:
                    semb_prior[j] = i[j]
                    times_cum[j] = a
                    times_max[j] = a*i[j]
                if semb > sembmax:
                    sembmax = semb
                    sembmaxX = x
                    sembmaxY = y
                    if times_min[j] == 0:
                        times_min[j] = a

            delta = orthodrome.distance_accurate50m_numpy(x, y, origin.lat,
                                                          origin.lon)
            azi = toAzimuth(float(Origin['lat']), float(Origin['lon']),
                            float(sembmaxX), float(sembmaxY))
            semb_prior = copy.copy(i)
            sembmaxvaluev[a] = sembmax
            sembmaxlatv[a] = sembmaxX
            sembmaxlonv[a] = sembmaxY
            fobjsembmax.write('%d %.3f %.3f %.30f %.30f %d %03f %f %03f\n'
                              % (a*step, sembmaxX, sembmaxY, sembmax, uncert,
                                 usedarrays, delta, float(azi), delta*119.19))
            fobj.close()

        fobjsembmax.close()

        fobj_cum = open(os.path.join(folder, 'semb_cum_%s_%s_boot%s_%s.ASC'
                                     % (switch, Origin['depth'], boot, phase)),
                        'w')
        for x, y, sembcums in zip(latv, lonv, semb_cum):
            fobj_cum.write('%.2f %.2f %.20f\n' % (x, y, sembcums))
        fobj_cum.close()

        fobj_timemax = open(os.path.join(folder, 'times_cum_%s_%s_boot%s_%s.ASC'
                                         % (switch, Origin['depth'], boot,
                                            phase)),
                            'w')
        for x, y, timemax in zip(latv, lonv, times_max):
            fobj_timemax.write('%.2f %.2f %.20f\n' % (x, y, timemax))
        fobj_timemax.close()

        fobj_timecum = open(os.path.join(folder, 'times_max_%s_%s_boot%s_%s.ASC'
                                         % (switch, Origin['depth'], boot,
                                            phase)),
                            'w')
        for x, y, timecum in zip(latv, lonv, times_cum):
            fobj_timecum.write('%.2f %.2f %.20f\n' % (x, y, timecum))
        fobj_timecum.close()

        fobj_timemin = open(os.path.join(folder, 'times_min_%s_%s_boot%s_%s.ASC'
                                         % (switch, Origin['depth'], boot,
                                            phase)),
                            'w')
        for x, y, timexy in zip(latv, lonv, times_min):
            fobj_timemin.write('%.2f %.2f %.20f\n' % (x, y, timexy))
        fobj_timemin.close()
    if cboot is None:
        fobjsembmax = open(os.path.join(folder, 'sembmax_%s_%s.txt'
                                        % (switch, phase)), 'w')
    else:
        fobjsembmax = open(os.path.join(folder, 'sembmax_%s_boot%s_%s.txt'
                                        % (switch, cboot, phase)), 'w')
    if temp_comb is not None:
        tmp_general = temp_comb
        tmp = tmp_general
        ishape = tmp_general[0]
        semb_cum = num.zeros(num.shape(ishape))
        times_cum = num.zeros(num.shape(ishape))
        times_min = num.zeros(num.shape(ishape))
        times_max = num.zeros(num.shape(ishape))
        semb_prior = num.zeros(num.shape(ishape))
    norm = num.max(num.max(tmp_general, axis=1))
    for a, i in enumerate(tmp_general):
        logger.info('timestep %d' % a)
        if cboot is None:
            fobj = open(os.path.join(folder, '%s-%s_%03d_%s.ASC'
                                     % (switch, Origin['depth'], a, phase)),
                        'w')
        else:
            fobj = open(os.path.join(folder, '%s-%s_boot%s_%03d_%s.ASC'
                                     % (switch, Origin['depth'], cboot, a,
                                        phase)),
                        'w')
        fobj.write('# %s , %s\n' % (d, rcs))
        fobj.write('# step %ds| ntimes %d| winlen: %ds\n'
                   % (step, ntimes, winlen))
        fobj.write('# \n')
        fobj.write('# southwestlat: %.2f dlat: %f nlat: %f \n'
                   % (Latul, gridspacing, dimX))
        fobj.write('# southwestlon: %.2f dlon: %f nlon: %f \n'
                   % (Lonul, gridspacing, dimY))
        fobj.write('# ddepth: 0 ndepth: 1 \n')

        sembmax = 0
        sembmaxX = 0
        sembmaxY = 0
        counter_time = 0
        uncert = num.std(i)
        semb_cum =+ i
        for j in range(migpoints):

            x = latv[j]
            y = lonv[j]

            if cfg.Bool('norm_all') is True:
                semb = i[j]/norm
            else:
                semb = i[j]
            fobj.write('%.2f %.2f %.20f\n' % (x, y, semb))
            if semb_prior[j] <= semb:
                semb_prior[j] = i[j]
                times_cum[j] = a
                times_max[j] = a*i[j]
            if semb > sembmax:
                sembmax = semb
                sembmaxX = x
                sembmaxY = y
                if times_min[j] == 0:
                    times_min[j] = a

        delta = orthodrome.distance_accurate50m_numpy(x, y, origin.lat,
                                                      origin.lon)
        azi = toAzimuth(float(Origin['lat']), float(Origin['lon']),
                        float(sembmaxX), float(sembmaxY))
        semb_prior = copy.copy(i)
        sembmaxvaluev[a] = sembmax
        sembmaxlatv[a] = sembmaxX
        sembmaxlonv[a] = sembmaxY
        fobjsembmax.write('%d %.3f %.3f %.30f %.30f %d %03f %f %03f\n'
                          % (a*step, sembmaxX, sembmaxY, sembmax, uncert,
                             usedarrays, delta, float(azi), delta*119.19))
        fobj.close()
    fobjsembmax.close()

    if cboot is None:
        fobj_cum = open(os.path.join(folder, 'semb_cum_%s_%s_%s.ASC'
                                     % (switch, Origin['depth'], phase)),
                        'w')
    else:
        fobj_cum = open(os.path.join(folder, 'semb_cum_%s_%s_boot%s_%s.ASC'
                                     % (switch, Origin['depth'], cboot, phase)),
                        'w')
    for x, y, sembcums in zip(latv, lonv, semb_cum):
        fobj_cum.write('%.2f %.2f %.20f\n' % (x, y, sembcums))
    fobj_cum.close()

    if cboot is None:
        fobj_timemax = open(os.path.join(folder, 'times_cum_%s_%s_%s.ASC'
                                         % (switch, Origin['depth'], phase)),
                            'w')
    else:
        fobj_timemax = open(os.path.join(folder, 'times_cum_%s_%s_boot%s_%s.ASC'
                                         % (switch, Origin['depth'], cboot,
                                            phase)),
                            'w')
    for x, y, timemax in zip(latv, lonv, times_max):
        fobj_timemax.write('%.2f %.2f %.20f\n' % (x, y, timemax))
    fobj_timemax.close()

    if cboot is None:
        fobj_timecum = open(os.path.join(folder, 'times_max_%s_%s_%s.ASC'
                                         % (switch, Origin['depth'], phase)),
                            'w')
    else:
        fobj_timecum = open(os.path.join(folder, 'times_max_%s_%s_boot%s_%s.ASC'
                                         % (switch, Origin['depth'], cboot,
                                            phase)),
                            'w')
    for x, y, timecum in zip(latv, lonv, times_cum):
        fobj_timecum.write('%.2f %.2f %.20f\n' % (x, y, timecum))
    fobj_timecum.close()

    if cboot is None:
        fobj_timemin = open(os.path.join(folder, 'times_min_%s_%s_%s.ASC'
                                         % (switch, Origin['depth'],
                                            phase)),
                            'w')
    else:
        fobj_timemin = open(os.path.join(folder, 'times_min_%s_%s_boot%s_%s.ASC'
                                         % (switch, Origin['depth'], cboot,
                                            phase)),
                            'w')
    for x, y, timexy in zip(latv, lonv, times_min):
        fobj_timemin.write('%.2f %.2f %.20f\n' % (x, y, timexy))
    fobj_timemin.close()

    trigger.writeSembMaxValue(sembmaxvaluev, sembmaxlatv, sembmaxlonv,
                              ntimes, Config, Folder)

    inspect_semb = cfg.Bool('inspect_semb')
    if inspect_semb is True and cboot is None:
        trigger.semblancestalta(sembmaxvaluev, sembmaxlatv, sembmaxlonv)
    return sembmaxvaluev, tmp


def collectSembweighted(SembList, Config, Origin, Folder, ntimes, arrays,
                        switch, weights):
    '''
    method to collect semblance matrizes from all processes and write them to file for each timestep
    '''
    Logfile.add('start collect in collectSemb')

    cfg= ConfigObj(dict=Config)
    origin = ConfigObj(dict=Origin)

    dimX = cfg.dimX()
    dimY = cfg.dimY()
    if switch == 0:
        winlen = cfg.winlen()
        step = cfg.step()
    if switch == 1:
        winlen = cfg.winlen_f2()
        step = cfg.step_f2()

    latv = []
    lonv = []

    gridspacing = cfg.Float('gridspacing')
    migpoints = dimX * dimY
    o_lat = origin.lat()
    o_lon = origin.lon()
    oLatul = 0
    oLonul = 0

    z=0

    for i in xrange(dimX):
         oLatul = o_lat -((dimX-1)/2) * gridspacing + i*gridspacing

         if z == 0 and i == 0 :
             Latul = oLatul
         o=0

         for j in xrange(dimY):
               oLonul = o_lon -((dimY-1)/2) * gridspacing + j*gridspacing

               if o==0 and j==0:
                    Lonul = oLonul

               latv.append(oLatul)
               lonv.append(oLonul)


    tmp=1
    weight_norm = num.sum(weights)
    for a, w in zip(SembList, weights):
        if num.mean(a)>0:
            tmp *= a*(w/weight_norm)

    sembmaxvaluev = num.ndarray(ntimes,dtype=float)
    sembmaxlatv = num.ndarray(ntimes,dtype=float)
    sembmaxlonv = num.ndarray(ntimes,dtype=float)

    rc= UTCDateTime(Origin['time'])
    rcs= '%s-%s-%s_%02d:%02d:%02d'%(rc.day,rc.month,rc.year, rc.hour,rc.minute,rc.second)
    d = rc.timestamp
    usedarrays = arrays


    folder = Folder['semb']
    fobjsembmax = open(os.path.join(folder,'sembmax_weighted_%s.txt' %(switch)),'w')

    for a, i in enumerate(tmp):
        logger.info('timestep %d' % a)


        fobj = open(os.path.join(folder,'%s-%s_%03d._weighted_semblance.ASC' %(switch,Origin['depth'],a)),'w')

        fobj.write('# %s , %s\n' %(d,rcs))
        fobj.write('# step %ds| ntimes %d| winlen: %ds\n' %(step,ntimes,winlen))
        fobj.write('# \n')
        fobj.write('# southwestlat: %.2f dlat: %f nlat: %f \n'%(Latul,gridspacing,dimX))
        fobj.write('# southwestlon: %.2f dlon: %f nlon: %f \n'%(Lonul,gridspacing,dimY))
        fobj.write('# ddepth: 0 ndepth: 1 \n')


        sembmax = 0
        sembmaxX = 0
        sembmaxY = 0

        origin = DataTypes.dictToLocation(Origin)
        uncert = num.std(i) #maybe not std?
        for j in range(migpoints):
            x= latv[j]
            y= lonv[j]
            semb = i[j]

            fobj.write('%.2f %.2f %.20f %.20f\n' %(x,y,semb))

            if  semb > sembmax:
                sembmax = semb# search for maximum and position of maximum on semblance grid for given time step
                sembmaxX = x
                sembmaxY = y

        delta = loc2degrees(Location(sembmaxX, sembmaxY), origin)
        azi = toAzimuth(float(Origin['lat']), float(Origin['lon']),float(sembmaxX), float(sembmaxY))

        sembmaxvaluev[a] = sembmax
        sembmaxlatv[a] = sembmaxX
        sembmaxlonv[a] = sembmaxY

        fobjsembmax.write('%d %.2f %.2f %.20f %.20f %d %03f %f %03f\n' %(a*step,sembmaxX,sembmaxY,sembmax,uncert,usedarrays,delta,float(azi),delta*119.19))
        fobj.close()


    fobjsembmax.close()

    trigger.writeSembMaxValue(sembmaxvaluev,sembmaxlatv,sembmaxlonv,ntimes,Config,Folder)
    trigger.semblancestalta(sembmaxvaluev,sembmaxlatv,sembmaxlonv)

def toMatrix(npVector, nColumns):

    t = npVector.tolist()[0]
    n = nColumns
    mat = []

    for i in range(len(t) / n):
        pos1 = i * n
        pos2 = pos1 + n
        slice = t[pos1:pos2]
        assert len(slice) == nColumns
        mat.append(slice)

    return mat


def doCalc(flag, Config, WaveformDict, FilterMetaData, Gmint, Gmaxt,
           TTTGridMap, Folder, Origin, ntimes, switch, ev, arrayfolder,
           syn_in, refshifts, phase, bs_weights=None):
    '''
    method for calculating semblance of one station array
    '''
    Logfile.add('PROCESS %d %s' % (flag,' Enters Semblance Calculation'))
    Logfile.add('MINT  : %f  MAXT: %f Traveltime' % (Gmint, Gmaxt))

    cfg = ConfigObj(dict=Config)
    cfg_f = FilterCfg(Config)

    timeev = util.str_to_time(ev.time)
    dimX = cfg.dimX()
    dimY = cfg.dimY()
    if switch == 0:
        winlen = cfg.Float('winlen')
        step = cfg.Float('step')
    if switch == 1:
        winlen = cfg.Float('winlen_f2')
        step = cfg.Float('step_f2')

    new_frequence = cfg.newFrequency()
    forerun = cfg.Float('forerun')
    duration = cfg.Float('duration')

    nostat = len(WaveformDict)
    traveltimes = {}
    recordstarttime = ''
    minSampleCount = 999999999

    if cfg.Float('forerun')>0:
        ntimes = float((cfg.Float('forerun') + cfg.Float('duration') ) / step )
    else:
        ntimes = float((cfg.Float('duration') ) / step )
    nsamp = float(winlen * new_frequence)
    nstep = float(step * new_frequence)

    obspy_compat.plant()

    ###########################################################################
    calcStreamMap = WaveformDict

    stations = []
    py_trs = []
    lats = []
    lons = []
    if cfg.Bool('synthetic_test') is False:
        if sys.version_info.major >= 3:
            for trace in sorted(calcStreamMap.keys()):
                py_tr = obspy_compat.to_pyrocko_trace(calcStreamMap[trace])
                py_trs.append(py_tr)
                for il in FilterMetaData:
                    if str(il) == str(trace):
                            szo = model.Station(lat=float(il.lat), lon=float(il.lon),
                                                station=il.sta, network=il.net,
                                                channels=py_tr.channel,
                                                elevation=il.ele, location=il.loc)
                            stations.append(szo)
                            lats.append(float(il.lat))
                            lons.append(float(il.lon))
        else:
            for trace in sorted(calcStreamMap.keys()):
                py_tr = obspy_compat.to_pyrocko_trace(calcStreamMap[trace])
                py_trs.append(py_tr)
                for il in FilterMetaData:
                    if str(il) == str(trace):
                            szo = model.Station(lat=float(il.lat), lon=float(il.lon),
                                                station=il.sta, network=il.net,
                                                channels=py_tr.channel,
                                                elevation=il.ele, location=il.loc)
                            stations.append(szo)
                            lats.append(float(il.lat))
                            lons.append(float(il.lon))
        array_center = [num.mean(lats), num.mean(lons)]

#==================================synthetic BeamForming======================

    if cfg.Bool('synthetic_test') is True:
        if sys.version_info.major >= 3:
            for trace in sorted(calcStreamMap.keys()):
                for il in FilterMetaData:
                    if str(il) == str(trace):
                            szo = model.Station(lat=float(il.lat), lon=float(il.lon),
                                                station=il.sta, network=il.net,
                                                channels='BHZ',
                                                elevation=il.ele, location=il.loc)
                            stations.append(szo)
                            lats.append(float(il.lat))
                            lons.append(float(il.lon))
        else:
            for trace in sorted(calcStreamMap.keys()):
                for il in FilterMetaData:
                    if str(il) == str(trace):
                            szo = model.Station(lat=float(il.lat), lon=float(il.lon),
                                                station=il.sta, network=il.net,
                                                channels='BHZ',
                                                elevation=il.ele, location=il.loc)
                            stations.append(szo)
                            lats.append(float(il.lat))
                            lons.append(float(il.lon))
        array_center = [num.mean(lats), num.mean(lons)]
        store_id = syn_in.store()
        engine = LocalEngine(store_superdirs=[syn_in.store_superdirs()])
        recordstarttimes = []
        targets = []
        sources = []
        for st in stations:
            target = Target(
                    lat=st.lat,
                    lon=st.lon,
                    store_id=store_id,
                    tmin=-1000,
                    tmax=1000,
                    codes=(st.network, st.station, st.location, 'BHZ'),
                    interpolation='multilinear',
                    quantity=cfg.quantity())
            targets.append(target)

        if syn_in.nsources() == 1:
            if syn_in.use_specific_stf() is True:
                stf = syn_in.stf()
                stf = exec(stf)
            else:
                stf = STF()
            if syn_in.source() == 'RectangularSource':
                    sources.append(RectangularSource(
                        lat=float(syn_in.lat_0()),
                        lon=float(syn_in.lon_0()),
                        east_shift=float(syn_in.east_shift_0())*1000.,
                        north_shift=float(syn_in.north_shift_0())*1000.,
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
                        time=util.str_to_time(syn_in.time_0())))
            if syn_in.source() == 'DCSource':
                    sources.append(DCSource(
                        lat=float(syn_in.lat_0()),
                        lon=float(syn_in.lon_0()),
                        east_shift=float(syn_in.east_shift_0())*1000.,
                        north_shift=float(syn_in.north_shift_0())*1000.,
                        depth=syn_in.depth_syn_0()*1000.,
                        strike=syn_in.strike_0(),
                        dip=syn_in.dip_0(),
                        rake=syn_in.rake_0(),
                        stf=stf,
                        time=util.str_to_time(syn_in.time_0()),
                        magnitude=syn_in.magnitude_0()))
            if syn_in.source() == 'MTSource':
                    sources.append(MTSource(
                        lat=float(syn_in.lat_0()),
                        lon=float(syn_in.lon_0()),
                        east_shift=float(syn_in.east_shift_0())*1000.,
                        north_shift=float(syn_in.north_shift_0())*1000.,
                        depth=syn_in.depth_syn_0()*1000.,
                        rmnn=syn_in.rmnn,
                        rmee=syn_in.rmee,
                        rmdd=syn_in.rmdd,
                        rmne=syn_in.rmne,
                        rmnd=syn_in.rmnd,
                        rmed=syn_in.rmed,
                        duation=syn_in.duration,
                        time=util.str_to_time(syn_in.time_0()),
                        magnitude=syn_in.magnitude_0()))
        else:
            for i in range(syn_in.nsources()):
                if syn_in.use_specific_stf() is True:
                    stf = syn_in.stf()
                    stf = exec(stf)

                else:
                    stf = STF()
                if syn_in.source() == 'RectangularSource':
                        sources.append(RectangularSource(
                            lat=float(syn_in.lat_1(i)),
                            lon=float(syn_in.lon_1(i)),
                            east_shift=float(syn_in.east_shift_1(i))*1000.,
                            north_shift=float(syn_in.north_shift_1(i))*1000.,
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
                            east_shift=float(syn_in.east_shift_1(i))*1000.,
                            north_shift=float(syn_in.north_shift_1(i))*1000.,
                            depth=syn_in.depth_syn_1(i)*1000.,
                            strike=syn_in.strike_1(i),
                            dip=syn_in.dip_1(i),
                            rake=syn_in.rake_1(i),
                            stf=stf,
                            time=util.str_to_time(syn_in.time_1(i)),
                            magnitude=syn_in.magnitude_1(i)))
        synthetic_traces = []
        for source in sources:
            response = engine.process(source, targets)
            synthetic_traces_source = response.pyrocko_traces()
            from pyrocko import trace as trld
            if not synthetic_traces:
                synthetic_traces = synthetic_traces_source
            else:
                for trsource, tr in zip(synthetic_traces_source,
                                        synthetic_traces):
                        tr.add(trsource)
            #debug
        #trld.snuffle(synthetic_traces)

        timeev = util.str_to_time(syn_in.time_0())
        if cfg.Bool('synthetic_test_add_noise') is True:
            from noise_addition import add_noise
            trs_orgs = []
            calcStreamMapsyn = calcStreamMap.copy()
            for tracex in calcStreamMapsyn.keys():
                    for trl in synthetic_traces:
                        if str(trl.name()[4:12])== str(tracex[4:]) or str(trl.name()[3:13])== str(tracex[3:]) or str(trl.name()[3:11])== str(tracex[3:]) or str(trl.name()[3:14])== str(tracex[3:]):
                            tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapsyn[tracex])
                            tr_org.downsample_to(2.0)
                            trs_orgs.append(tr_org)
            store_id = syn_in.store()
            engine = LocalEngine(store_superdirs=[syn_in.store_superdirs()])
            synthetic_traces = add_noise(trs_orgs, engine, source.pyrocko_event(),
                                         stations,
                                         store_id, phase_def=phase)
        trs_org = []
        trs_orgs = []
        from pyrocko import trace
        fobj = os.path.join(arrayfolder, 'shift.dat')
        calcStreamMapsyn = {}
        if cfg.Bool('synthetic_test_pertub_arrivals') is True:
            shift_max = cfg.Float('shift_max')
            for trl in synthetic_traces:
                shift = num.random.uniform(-shift_max,shift_max)
                trl.shift(shift)

        if sys.version_info.major >= 3:
            for tracex, trl in zip(sorted(calcStreamMap.keys()), synthetic_traces):
                        if cfg.Bool('dynamic_filter') is False:
                            if switch == 0:
                                trl.bandpass(4,cfg_f.flo(), cfg_f.fhi())
                            elif switch == 1:
                                trl.bandpass(4,cfg_f.flo2(), cfg_f.fhi2())
                        synthetic_obs_tr = obspy_compat.to_obspy_trace(trl)
                        calcStreamMap[tracex] = synthetic_obs_tr
        else:
            for tracex, trl in zip(calcStreamMap.keys(), synthetic_traces):
                            if cfg.Bool('dynamic_filter') is False:
                                if switch == 0:
                                    trl.bandpass(4,cfg_f.flo(), cfg_f.fhi())
                                elif switch == 1:
                                    trl.bandpass(4,cfg_f.flo2(), cfg_f.fhi2())
                            synthetic_obs_tr = obspy_compat.to_obspy_trace(trl)
                            calcStreamMap[tracex] = synthetic_obs_tr

    if cfg.Bool('shift_by_phase_pws') is True:
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

    if cfg.Bool('shift_by_phase_cc') is True:
        from stacking import align_traces
        calcStreamMapshifted= calcStreamMap.copy()
        list_tr = []
        for trace in calcStreamMapshifted.keys():
            tr_org = calcStreamMapshifted[trace]
            list_tr.append(tr_org)
        shifts, ccs = align_traces(list_tr, 10, master=False)
        for shift in shifts:
            for trace in calcStreamMapshifted.keys():
                    tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapshifted[trace])
                    tr_org.shift(shift)
                    shifted = obspy_compat.to_obspy_trace(tr_org)
                    calcStreamMapshifted[trace] = shifted
        calcStreamMap = calcStreamMapshifted

    if cfg.Bool('shift_by_phase_onset') is True:
        pjoin = os.path.join
        timeev = util.str_to_time(ev.time)
        trs_orgs = []
        calcStreamMapshifted = calcStreamMap.copy()
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
        for tracex in calcStreamMapshifted.keys():
                for trl in shifted_traces:
                    if str(trl.name()[4:12]) == str(tracex[4:]) or str(trl.name()[3:13])== str(tracex[3:]) or str(trl.name()[3:11])== str(tracex[3:]) or str(trl.name()[3:14])== str(tracex[3:]):
                        mod = trl
                        recordstarttime = calcStreamMapshifted[tracex].stats.starttime.timestamp
                        recordendtime = calcStreamMapshifted[tracex].stats.endtime.timestamp
                        tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapshifted[tracex])
                        tr_org_add = mod.chop(recordstarttime, recordendtime, inplace=False)
                        shifted_obs_tr = obspy_compat.to_obspy_trace(tr_org_add)
                        calcStreamMapshifted[tracex] = shifted_obs_tr
        calcStreamMap = calcStreamMapshifted

    weight = 1.
    if cfg.Bool('weight_by_noise') is True:
        from noise_analyser import analyse
        pjoin = os.path.join
        timeev = util.str_to_time(ev.time)
        trs_orgs = []
        calcStreamMapshifted = calcStreamMap.copy()
        for trace in calcStreamMapshifted.keys():
                tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapshifted[trace])
                trs_orgs.append(tr_org)

        timing = CakeTiming(
           phase_selection='first(p|P|PP|P(cmb)P(icb)P(icb)p(cmb)p)-20',
           fallback_time=100.)
        traces = trs_orgs
        event = model.Event(lat=float(ev.lat), lon=float(ev.lon),
                            depth=ev.depth*1000., time=timeev)
        directory = arrayfolder
        bf = BeamForming(stations, traces, normalize=True)
        shifted_traces = bf.process(event=event,
                                    timing=timing,
                                    fn_dump_center=pjoin(directory,
                                                         'array_center.pf'),
                                    fn_beam=pjoin(directory, 'beam.mseed'))
        i = 0
        store_id = syn_in.store()
        engine = LocalEngine(store_superdirs=[syn_in.store_superdirs()])
        weight = analyse(shifted_traces, engine, event, stations,
                         100., store_id, nwindows=1,
                         check_events=True, phase_def=phase)

    try:
        if cfg.Bool('futterman_attenuation') is True and phase is 'S':
            trs_orgs = []
            for trace in calcStreamMap.keys():
                    tr_org = obspy_compat.to_pyrocko_trace(calcStreamMap[trace])
                    mod = cake.load_model(crust2_profile=(ev.lat, ev.lon))
                    timing = CakeTiming(
                       phase_selection='first(S)',
                       fallback_time=100.)
                    dist = ortho.distance_accurate50m(ev, event1)
                    ray = timing.t(mod,(ev.depth, dist), get_ray=True)
                    zx, xx, tx = ray.zxt_path_subdivided()
                    qs_int = 0.
                    vs_int = 0.
                    vp_int = 0.
                    qp_int = 0.

                    zx = zx[0]
                    for z in zx:
                    	mat = mod.material(z)
                    	qs_int = qs_int + mat.qs
                    	vs_int = vs_int + mat.vs
                    	qp_int = qp_int + mat.qp
                    	vp_int = vp_int + mat.vp
                    #dist =ray.x*(d2r*earthradius/km)
                    #dtstar = dist/(qs_int*vs_int) #direct dstar
                    L = 200.e3
                    dtstar = L*(1./vp_int/qp_int - 1./vs_int/qs_int) # differential tstar measurement
                    npts = len(tr_org.ydata)
                    idat = tr_org.ydata
                    ffs = num.fft.rfft(idat)

                    ff, IDAT = spectrum(tr_org.ydata, tr_org.deltat)

                    w = 2*num.pi*ff
                    wabs = abs(w)

                    w0 = 2*num.pi
                    Aw = num.exp(-0.5*wabs*dtstar)
                    phiw = (1/num.pi)*dtstar*num.log(2*num.pi*num.exp(num.pi)/wabs)
                    # phiw = 0.5 * (wabs/w0)**(-1) * dtstar * 1./num.tan(1*num.pi/2)
                    phiw = num.nan_to_num(phiw)
                    Dwt = Aw * num.exp(-1j*w*phiw)
                    qdat = num.real(num.fft.ifft((IDAT*Dwt)))
                    tr_org.ydata = qdat
                    trs_orgs.append(tr_org)
            for tracex in calcStreamMap.keys():
                    for trl in trs_orgs:
                        obs_tr = obspy_compat.to_obspy_trace(trl)
                        calcStreamMap[tracex] = obs_tr
    except:
        pass

    if cfg.Bool('array_response') is True:
        from obspy.signal import array_analysis
        from obspy.core import stream
        ntimesr = float((forerun + duration)/step)
        nsampr = float(winlen)
        nstepr = float(step)
        sll_x=-3.0
        slm_x=3.0
        sll_y=-3.0
        slm_y=3.0
        sl_s=0.03,
        # sliding window properties

        # frequency properties
        frqlow=1.0,
        frqhigh=8.0
        prewhiten=0
        # restrict output
        semb_thres=-1e9
        vel_thres=-1e9
        stime=stime
        etime=etime
        stream_arr = stream.Stream()
        for trace in calcStreamMapshifted.keys():
            stream_arr.append(calcStreamMapshifted[trace])
        results = array_analysis.array_processing(stream_arr, nsamp, nstep,
                                                  sll_x, slm_x, sll_y, slm_y,
                                                   sl_s, semb_thres, vel_thres,
                                                   frqlow, frqhigh, stime,
                                                   etime, prewhiten)
        timestemp = results[0]
        relative_relpow = results[1]
        absolute_relpow = results[2]
    if sys.version_info.major >= 3:
        for trace in sorted(calcStreamMap.keys()):
            recordstarttime = calcStreamMap[trace].stats.starttime
            d = calcStreamMap[trace].stats.starttime
            d = d.timestamp

            if calcStreamMap[trace].stats.npts < minSampleCount:
                minSampleCount = calcStreamMap[trace].stats.npts
    else:
        for trace in calcStreamMap.keys():
            recordstarttime = calcStreamMap[trace].stats.starttime
            d = calcStreamMap[trace].stats.starttime
            d = d.timestamp

            if calcStreamMap[trace].stats.npts < minSampleCount:
                minSampleCount = calcStreamMap[trace].stats.npts

    traces = num.ndarray(shape=(len(calcStreamMap), minSampleCount),
                         dtype=float)
    traveltime = num.ndarray(shape=(len(calcStreamMap), dimX*dimY),
                             dtype=float)

    latv = num.ndarray(dimX*dimY, dtype=float)
    lonv = num.ndarray(dimX*dimY, dtype=float)

    c = 0
    streamCounter = 0
    if sys.version_info.major >= 3:
        for key in sorted(calcStreamMap.keys()):
            streamID = key
            c2 = 0

            for o in calcStreamMap[key]:
                if c2 < minSampleCount:
                    traces[c][c2] = o

                    c2 += 1

            for key in sorted(TTTGridMap.keys()):

                if streamID == key:
                    traveltimes[streamCounter] = TTTGridMap[key]
                else:
                    "NEIN", streamID, key

            if streamCounter not in traveltimes:
                continue

            g = traveltimes[streamCounter]
            dimZ = g.dimZ
            mint = g.mint
            gridElem = g.GridArray

            for x in range(dimX):
                for y in range(dimY):
                    elem = gridElem[x, y]

                    traveltime[c][x * dimY + y] = elem.tt
                    latv[x * dimY + y] = elem.lat
                    lonv[x * dimY + y] = elem.lon
            c += 1
            streamCounter += 1

    else:
        for key in calcStreamMap.keys():
            streamID = key
            c2 = 0

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
               continue

            g = traveltimes[streamCounter]
            dimZ = g.dimZ
            mint = g.mint
            gridElem = g.GridArray

            for x in range(dimX):
                for y in range(dimY):
                    elem = gridElem[x, y]

                    traveltime [c][x * dimY + y] = elem.tt
                    latv [x * dimY + y] = elem.lat
                    lonv [x * dimY + y] = elem.lon

            c += 1
            streamCounter += 1


    ################ CALCULATE PARAMETER FOR SEMBLANCE CALCULATION ########
    nsamp = winlen * new_frequence

    nstep = step*new_frequence
    migpoints = dimX * dimY

    dimZ = 0
    maxp = int(Config['ncore'])

    Logfile.add('PROCESS %d  NTIMES: %d' %(flag,ntimes))

    if False:
        print('nostat ',nostat,type(nostat))
        print('nsamp ',nsamp,type(nsamp))
        print('ntimes ',ntimes,type(ntimes))
        print('nstep ',nstep,type(nstep))
        print('dimX ',dimX,type(dimX))
        print('dimY ',dimY,type(dimY))
        print('mint ',Gmint,type(mint))
        print('new_freq ',new_frequence,type(new_frequence))
        print('minSampleCount ',minSampleCount,type(minSampleCount))
        print('latv ',latv,type(latv))
        print('traces',traces,type(traces))

#===================compressed sensing=================================
    try:
        cs = cfg.cs()
    except:
        cs = 0
    if cs == 1:
        csmaxvaluev = num.ndarray(ntimes,dtype=float)
        csmaxlatv = num.ndarray(ntimes,dtype=float)
        csmaxlonv = num.ndarray(ntimes,dtype=float)
        folder = Folder['semb']
        fobjcsmax = open(os.path.join(folder,'csmax_%s.txt' %(switch)),'w')
        traveltimes = traveltime.reshape(1,nostat*dimX*dimY)
        traveltime2 = toMatrix(traveltimes, dimX * dimY)  # for relstart
        traveltime = traveltime.reshape(dimX*dimY,nostat)
        import matplotlib as mpl
        import scipy.optimize as spopt
        import scipy.fftpack as spfft
        import scipy.ndimage as spimg
        import cvxpy as cvx
        import matplotlib.pyplot as plt
        A = spfft.idct(traveltime, norm='ortho', axis=0)
        n =(nostat*dimX*dimY)
        vx = cvx.Variable(dimX*dimY)
        res = cvx.Variable(1)
        objective = cvx.Minimize(cvx.norm(res, 1))
        back2 = num.zeros([dimX,dimY])
        l = int(nsamp)
        fobj = open(os.path.join(folder,'%s-%s_%03d.cs' %(switch,Origin['depth'],l)),'w')
        for i in range(ntimes):
            ydata = []
            try:
                for tr in traces:
                        relstart = int((dimX*dimY - mint) * new_frequence + 0.5) + i * nstep
                        tr=spfft.idct(tr[relstart+i:relstart+i+dimX*dimY], norm='ortho', axis=0)

                        ydata.append(tr)
                        ydata = num.asarray(ydata)
                        ydata = ydata.reshape(dimX*dimY,nostat)

                        constraints = [res == cvx.sum_entries( 0+ num.sum([ydata[:,x]-A[:,x]*vx  for x in range(nostat) ]) ) ]

                        prob = cvx.Problem(objective, constraints)
                        result = prob.solve(verbose=False, max_iters=200)

                        x = num.array(vx.value)
                        x = num.squeeze(x)
                        back1 = x.reshape(dimX,dimY)
                        sig = spfft.idct(x, norm='ortho', axis=0)
                        back2 = back2 + back1
                        xs = num.array(res.value)
                        xs = num.squeeze(xs)
                        max_cs = num.max(back1)
                        idx = num.where(back1==back1.max())
                        csmaxvaluev[i] = max_cs
                        csmaxlatv[i] = latv[idx[0]]
                        csmaxlonv[i] = lonv[idx[1]]
                        fobj.write('%.5f %.5f %.20f\n' %(latv[idx[0]],lonv[idx[1]],max_cs))
                        fobjcsmax.write('%.5f %.5f %.20f\n' %(latv[idx[0]],lonv[idx[1]],max_cs))
                fobj.close()
                fobjcsmax.close()

            except:
                pass


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
        if cfg.Float('forerun') > 0:
            ntimes = int((cfg.Float('forerun') + cfg.Float('duration'))/step)
        else:
            ntimes = int((cfg.Float('duration')) / step)
        nsamp = float(winlen)
        nstep = float(step)
        Gmint = cfg.Int('forerun')

        k = semblance(maxp, nostat, nsamp, ntimes, nstep, dimX, dimY, Gmint,
                      new_frequence, minSampleCount, latv, lonv, traveltimes,
                      traces, calcStreamMap, timeev, Config, Origin, refshifts,
                      bs_weights=bs_weights)
        print("--- %s seconds ---" % (time.time() - start_time))

    t2 = time.time()

    Logfile.add('%s took %0.3f s' % ('CALC:', (t2-t1)))

    partSemb = k
    partSemb = partSemb.reshape(ntimes, migpoints)

    return partSemb, weight, array_center
