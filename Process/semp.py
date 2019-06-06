import os
import sys
import time
import numpy as num
from threading import Thread
from obspy.signal.headers import clibsignal
from obspy import Stream, Trace
# add local directories to import path
sys.path.append('../Common/')
from ConfigFile import ConfigObj, OriginCfg, SynthCfg, FilterCfg

import numpy as np

import Logfile
import Basic
import ctypes as C
from pyrocko import trace as trld
from pyrocko.marker import PhaseMarker

trace_txt  = 'trace.txt'
travel_txt = 'travel.txt'
latv_txt   = 'latv.txt'
lonv_txt   = 'lonv.txt'
semb_txt   = 'semb.txt'

def normalize(lst):
    s = sum(lst)
    return map(lambda x: float(x)/s, lst)


def hilbert(x, N=None):
    '''
    Return the hilbert transform of x of length N.
    (from scipy.signal, but changed to use fft and ifft from numpy.fft)
    '''
    x = num.asarray(x)
    if N is None:
        N = len(x)
    if N <= 0:
        raise ValueError("N must be positive.")
        x = num.real(x)
    Xf = num.fft.fft(x, N, axis=0)
    h = num.zeros(N)
    if N % 2 == 0:
        h[0] = h[N//2] = 1
        h[1:N//2] = 2
    else:
        h[0] = 1
        h[1:(N+1)//2] = 2

    if len(x.shape) > 1:
        h = h[:, num.newaxis]
    x = num.fft.ifft(Xf*h)
    return x


def xcorr(tr1, tr2, shift_len, full_xcorr=False):

    from obspy.core.util.deprecation_helpers import ObsPyDeprecationWarning
    if min(len(tr1), len(tr2)) - 2 * shift_len <= 0:
        msg = "shift_len too large. The underlying C code would silently " + \
              "use shift_len/2 which we want to avoid."
        raise ValueError(msg)
    tr1 = np.ascontiguousarray(tr1, np.float32)
    tr2 = np.ascontiguousarray(tr2, np.float32)
    corp = np.empty(2 * shift_len + 1, dtype=np.float64, order='C')

    shift = C.c_int()
    coe_p = C.c_double()

    res = clibsignal.X_corr(tr1, tr2, corp, shift_len, len(tr1), len(tr2),
                            C.byref(shift), C.byref(coe_p))
    if res:
        raise MemoryError

    if full_xcorr:
        return shift.value, coe_p.value, corp
    else:
        return shift.value, coe_p.value



class MyThread(Thread):

    def __init__(self, nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount):
        Thread.__init__(self)

        self.nostat= nostat
        self.nsamp = nsamp
        self.i= i
        self.nstep = nstep
        self.dimX  = dimX
        self.dimY  = dimY
        self.mint  = mint
        self.new_freq  = new_freq
        self.minSampleCount = minSampleCount


def toMatrix(npVector, nColumns):

    t   = npVector.tolist()[0]
    n   = nColumns
    mat = []

    for i in range(int(len(t) / n)) :
       pos1  = i * n
       pos2  = pos1 + n
       slice = t [pos1:pos2]
       assert len(slice) == nColumns
       mat.append(slice)

    return mat


def semblance(ncpus, nostat, nsamp, ntimes, nstep, dimX, dimY, mint,
              new_frequence, minSampleCount, latv_1, lonv_1, traveltime_1,
              trace_1, calcStreamMap, time, Config, Origin, refshifts,
              bs_weights=None):

        cfg = ConfigObj(dict=Config)
        origin = OriginCfg(Origin)
        cfg_f = FilterCfg(Config)

        if cfg.Bool('dynamic_filter') is False:

           return semblance_py(ncpus, nostat, nsamp, ntimes, nstep, dimX, dimY,
                               mint, new_frequence, minSampleCount, latv_1,
                               lonv_1, traveltime_1, trace_1, calcStreamMap,
                               time, cfg, refshifts, bs_weights=bs_weights)
        else:
           return semblance_py_dynamic_cf(ncpus, nostat, nsamp, ntimes, nstep,
                                          dimX, dimY, mint, new_frequence,
                                          minSampleCount, latv_1, lonv_1,
                                          traveltime_1, trace_1, calcStreamMap,
                                          time, origin, cfg_f)


def t2ind_fast(t, tdelta, snap=round):
    return int(int((t/tdelta)*(10**0))/(10.**0))

def t2ind(t, tdelta, snap=round):
    return int(snap(t/tdelta))


def semblance_py_dynamic_cf(ncpus, nostat, nsamp, ntimes, nstep, dimX, dimY,
                            mint, new_frequence, minSampleCount, latv_1, lonv_1,
                            traveltime_1, trace_1, calcStreamMap, time,
                            origin, FilterCfg):

    from pyrocko import obspy_compat
    obspy_compat.plant()
    trs_orgs  = []
    for tr in calcStreamMap:
        tr_org = obspy_compat.to_pyrocko_trace(calcStreamMap[tr])
        tr_org.ydata = tr_org.ydata / np.sqrt(np.mean(np.square(tr_org.ydata)))
        trs_orgs.append(tr_org)
    trace  = toMatrix(trace_1, minSampleCount)
    traveltime = []
    traveltime = toMatrix(traveltime_1, dimX * dimY)

    latv   = latv_1.tolist()
    lonv   = lonv_1.tolist()

    '''
    Basic.writeMatrix(trace_txt,  trace, nostat, minSampleCount, '%e')
    Basic.writeMatrix(travel_txt, traveltime, nostat, dimX * dimY, '%e')
    Basic.writeVector(latv_txt,   latv, '%e')
    Basic.writeVector(lonv_txt,   lonv, '%e')
    '''
    snap=(round, round)
    backSemb = np.ndarray(shape=(ntimes, dimX*dimY), dtype=float)
    for i in range(ntimes):
        #  loop over grid points
        sembmax = 0; sembmaxX = 0; sembmaxY = 0

        for j in range(dimX * dimY):
            semb = 0; nomin = 0; denom = 0
            sums_cc = 0
            sums = 0
            shifted = []
            relstart = []
            relstarts = nostat
            cc_data = []
            tt = []

            for k in range(nostat):
                relstart = traveltime[k][j]
                tr = trs_orgs[k]
                tmin = time+relstart+(i*nstep)-mint
                tmax = time+relstart+(i*nstep)-mint+nsamp
                try:
                    ibeg = max(0, t2ind_fast(tmin-tr.tmin, tr.deltat, snap[0]))
                    iend = min(
                        tr.data_len(),
                        t2ind_fast(tmax-tr.tmin, tr.deltat, snap[1]))
                except:
                    print('Loaded traveltime grid wrong!')

                data = tr.ydata[ibeg:iend]
                try:
                    sums += num.gradient(abs(data))
                except:
                    pass
                relstarts -=(relstart)

            sum = abs(num.sum(((sums))))
            denom = sum**2
            nomin = sum
            semb = sum

            backSemb[i][j] = sum
            if semb > sembmax :

               sembmax  = semb
               sembmaxX = latv[j]
               sembmaxY = lonv[j]

        Logfile.add('max semblance: ' + str(sembmax) + ' at lat/lon: ' +
                     str(sembmaxX)+','+ str(sembmaxY))

    backSemb = backSemb/num.max(num.max(backSemb))
    return abs(backSemb)


def semblance_py(ncpus, nostat, nsamp, ntimes, nstep, dimX, dimY, mint,
                 new_frequence, minSampleCount, latv_1, lonv_1, traveltime_1,
                 trace_1, calcStreamMap, time, cfg, refshifts,
                 bs_weights=None):
    from pyrocko import obspy_compat
    obspy_compat.plant()
    trs_orgs = []
    snap = (round, round)

    for tr in sorted(calcStreamMap):
        tr_org = obspy_compat.to_pyrocko_trace(calcStreamMap[tr])
        tr_org.ydata = tr_org.ydata / np.sqrt(np.mean(np.square(tr_org.ydata)))
        if cfg.Bool('combine_all') is True:
            # some trickery to make all waveforms have same polarity, while still
            # considering constructive/destructive interferences. This is needed
            # when combing all waveforms/arrays from the world at once(only then)
            # for a single array with polarity issues we recommend fixing polarity.
            # advantage of the following is that nothing needs to be known about the
            # mechanism.
            tr_org.ydata = abs(tr_org.ydata)
            tr_org.ydata = num.diff(tr_org.ydata)
        trs_orgs.append(tr_org)

    traveltime = []
    traveltime = toMatrix(traveltime_1, dimX * dimY)
    latv = latv_1.tolist()
    lonv = lonv_1.tolist()
    '''
    Basic.writeMatrix(trace_txt,  trace, nostat, minSampleCount, '%e')
    Basic.writeMatrix(travel_txt, traveltime, nostat, dimX * dimY, '%e')
    Basic.writeVector(latv_txt,   latv, '%e')
    Basic.writeVector(lonv_txt,   lonv, '%e')
    '''
    if nsamp == 0:
        nsamp = 1
    backSemb = np.ndarray(shape=(ntimes, dimX*dimY), dtype=float)
    for i in range(ntimes):
        sembmax = 0
        sembmaxX = 0
        sembmaxY = 0
        for j in range(dimX * dimY):
            semb = 0
            nomin = 0
            denom = 0
            if cfg.Bool('combine_all') is True:
                sums = 1
            else:
                sums = 0
            relstart = []
            relstarts = nostat

            for k in range(nostat):
                relstart = traveltime[k][j]
                tr = trs_orgs[k]
                tmin = time+relstart+(i*nstep)-mint-refshifts[k]
                tmax = time+relstart+(i*nstep)-mint+nsamp-refshifts[k]
                try:
                    ibeg = max(0, t2ind_fast(tmin-tr.tmin, tr.deltat, snap[0]))
                    iend = min(
                        tr.data_len(),
                        t2ind_fast(tmax-tr.tmin, tr.deltat, snap[1]))
                except Exception:
                    print('Loaded traveltime grid wrong!')

                data = tr.ydata[ibeg:iend]

                try:
                    if cfg.Bool('combine_all') is True:
                        if cfg.Bool('bootstrap_array_weights') is True:
                            sums *= (data)*bs_weights[k]
                        else:
                            sums *= (data)

                    else:
                        sums += (data)
                except Exception:
                    pass
                relstarts -= relstart

            sum = abs(num.sum(((sums))))
            semb = sum

            backSemb[i][j] = sum
            if semb > sembmax:
                sembmax  = semb   # search for maximum and position of maximum on semblance
                                 # grid for given time step
                sembmaxX = latv[j]
                sembmaxY = lonv[j]
        Logfile.add('max semblance: ' + str(sembmax) + ' at lat/lon: ' +
                     str(sembmaxX)+','+ str(sembmaxY))

    backSemb = backSemb
    return abs(backSemb)


def execsemblance(nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount) :

    f = [nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount]
    args  = Basic.floatToString(f, delim= ',')
    prog  = sys.executable + ' ' + __file__
    cmd   = prog  + ' ' + args

    Logfile.add('--------------------------------------------', cmd)
    result = Basic.systemCmd(cmd)
    Logfile.addLines(result)
    backSemb = Basic.readVector(semb_txt)
    return backSemb

def execsemblance2() :

    for i in range(len(sys.argv)):
        print(sys.argv[i])

    params = Basic.stringToFloat(sys.argv[1])
    [nostat, nsamp, i, nstep, dimX, dimY, mint, new_freq, minSampleCount] = params
    backSemb = startsemblance(int(nostat), int(nsamp), int(i), int(nstep),
                              int(dimX), int(dimY),  mint, new_freq,
                              int(minSampleCount), False)

    Basic.writeVector(semb_txt, backSemb)


def startsemblance(nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount, isParent = True):

    backSemb = []

    if isParent:
        backSemb = execsemblance(nostat, nsamp, i, nstep, dimX, dimY, mint,
                                 new_freq, minSampleCount)

    else:
       trace = Basic.readMatrix(trace_txt,  nostat, minSampleCount, '%e')
       traveltime = Basic.readMatrix(travel_txt, nostat, dimX * dimY, '%e')
       latv = Basic.readVector(latv_txt, '%e')
       lonv = Basic.readVector(lonv_txt, '%e')

    for j in range(dimX * dimY):
      semb  = 0
      nomin = 0
      denom = 0

      for l in range(int(nsamp)):
         sum = 0

         for k in range(nostat):
            relstart_samples = int((traveltime[k][j] - mint) * new_freq + 0.5) + i * nstep

            val   =  trace[k][relstart_samples+l]
            sum   += val
            denom +=(val * val)

         nomin += sum * sum;
         semb  = nomin /(float(nostat) * denom)

      backSemb.append(semb)


    return backSemb

# -------------------------------------------------------------------------------------------------

if __name__ == "__main__":
   execsemblance2()
