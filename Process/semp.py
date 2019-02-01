import os
import sys
import time
import numpy as num
from threading import Thread
from obspy.signal.headers import clibsignal
from obspy import Stream, Trace
# add local directories to import path
sys.path.append ('../Common/')

import numpy as np          #    Import from NumPy

import Logfile              #    Import from common
import Basic
import ctypes as C

trace_txt  = 'trace.txt'
travel_txt = 'travel.txt'
latv_txt   = 'latv.txt'
lonv_txt   = 'lonv.txt'
semb_txt   = 'semb.txt'

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

def startC_Code (nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_freq, minSampleCount) :

        f= [nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_freq, minSampleCount]
        args = Basic.floatToString (f, delim= ',')
        path =  os.path.dirname (__file__)

        prog = os.path.join (path,'OTest')

        os.system (prog  + ' ' + args)

# -------------------------------------------------------------------------------------------------

class MyThread (Thread):

    def __init__ (self, nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount):
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

    def run(self):

        return startC_Code (self.nostat, self.nsamp, self.i, self.nstep, self.dimX,self.dimY,
                            self.mint, self.new_freq, self.minSampleCount)

# -------------------------------------------------------------------------------------------------

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


def otest (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_frequence, minSampleCount,
               latv_1, lonv_1, traveltime_1, trace_1, calcStreamMap, time) :
   #if USE_C_CODE  :
    USE_C_CODE = False
    if USE_C_CODE == True :
       return otestSeriell (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_frequence,
                            minSampleCount, latv_1, lonv_1, traveltime_1, trace_1)
    else :
       return otest_py   (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_frequence,
                            minSampleCount, latv_1, lonv_1, traveltime_1, trace_1, calcStreamMap, time)

def t2ind_fast(t, tdelta, snap=round):
    return int(int((t/tdelta)*(10**0))/(10.**0))

def t2ind(t, tdelta, snap=round):
    return int(snap(t/tdelta))

def otest_py(ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_frequence, minSampleCount,
               latv_1, lonv_1, traveltime_1, trace_1, calcStreamMap, time) :
    from pyrocko import obspy_compat
    obspy_compat.plant()
    trs_orgs  = []
    for tr in calcStreamMap:
        trs_orgs.append(tr_org)
        tr_org = obspy_compat.to_pyrocko_trace(calcStreamMap[tr])
    trace  = toMatrix (trace_1, minSampleCount)
    traveltime = []
    traveltime = toMatrix (traveltime_1, dimX * dimY)

    latv   = latv_1.tolist()
    lonv   = lonv_1.tolist()

    '''
    Basic.writeMatrix (trace_txt,  trace, nostat, minSampleCount, '%e')
    Basic.writeMatrix (travel_txt, traveltime, nostat, dimX * dimY, '%e')
    Basic.writeVector (latv_txt,   latv, '%e')
    Basic.writeVector (lonv_txt,   lonv, '%e')
    '''
    snap= (round, round)

    backSemb = np.ndarray (shape=(ntimes, dimX*dimY), dtype=float)
    for i in range (ntimes) :
        #  loop over grid points
        sembmax = 0; sembmaxX = 0; sembmaxY = 0

        for j in range (dimX * dimY):
            semb = 0; nomin = 0; denom = 0
            sums_cc = 0
            sums = 0
            shifted = []
            relstart = []
            relstarts = nostat
            cc_data = []
            for k in range (nostat):
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
                    sums += (data) ##put gradient on top
                except:
                    pass
                relstarts -= (relstart)

            #for dat in cc_data:
        #        sums_cc +=xcorr(cc_data[0],dat,0)[1]
            sum = abs(num.sum(((sums))))
        #    sum = sums_cc
            denom = sum**2
            nomin = sum
            semb = sum

            backSemb[i][j] = sum
            if semb > sembmax :

               sembmax  = semb   # search for maximum and position of maximum on semblance
                                 # grid for given time step
               sembmaxX = latv[j]
               sembmaxY = lonv[j]

        Logfile.add ('max semblance: ' + str(sembmax) + ' at lat/lon: ' +
                     str(sembmaxX)+','+ str (sembmaxY))

    backSemb = backSemb/num.max(num.max(backSemb))

    return abs(backSemb)

# -------------------------------------------------------------------------------------------------

def otestSeriell (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_freq, minSampleCount,
                  latv_1, lonv_1, traveltime_1, trace_1) :

    trace  = toMatrix (trace_1, minSampleCount)
    traveltime = toMatrix (traveltime_1, dimX * dimY)

    Basic.writeMatrix (trace_txt,  trace, nostat, minSampleCount)
    Basic.writeMatrix (travel_txt, traveltime, nostat, dimX * dimY)
    Basic.writeVector (latv_txt,   latv_1.tolist())
    Basic.writeVector (lonv_txt,   lonv_1.tolist())
    '''
    Basic.writeMatrix (trace_txt,  trace, nostat, minSampleCount, '%e')
    Basic.writeMatrix (travel_txt, traveltime, nostat, dimX * dimY, '%e')
    Basic.writeVector (latv_txt,   latv_1.tolist(), '%e')
    Basic.writeVector (lonv_txt,   lonv_1.tolist(), '%e')
    '''

    startC_Code (nostat, int (nsamp), ntimes, nstep, dimX,dimY, mint, new_freq, minSampleCount)

    result   = Basic.readMatrix (semb_txt, ntimes, dimX*dimY)
    backSemb = np.ndarray (shape=(ntimes, dimX*dimY), dtype=float)

    for i in range (ntimes) :
        for j in range (dimX * dimY) :
            backSemb [i][j] = result [i][j]
    return backSemb

# -------------------------------------------------------------------------------------------------

def otestPar (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_freq, minSampleCount,
              latv_1, lonv_1, traveltime_1, trace_1) :

    trace  = toMatrix (trace_1, minSampleCount)
    traveltime = toMatrix (traveltime_1, dimX * dimY)
    latv   = latv_1.tolist()
    lonv   = lonv_1.tolist()

    Basic.writeMatrix (trace_txt,  trace, nostat, minSampleCount, '%e')
    Basic.writeMatrix (travel_txt, traveltime, nostat, dimX * dimY, '%e')
    Basic.writeVector (latv_txt,   latv, '%e')
    Basic.writeVector (lonv_txt,   lonv, '%e')

    backSemb = np.ndarray (shape=(ntimes, dimX*dimY), dtype=float)
    threads  = []

    for i in range(ntimes) :
        t = MyThread (nostat, int(nsamp), i, nstep, dimX,dimY, mint, new_freq, minSampleCount)
        t.run()
        threads.append (t)

    for t in threads :
        t.join()


    for i in range (ntimes) :
       #  loop over grid points
       #
       sembmax = 0
       sembmaxX= 0
       sembmaxY= 0
       backSemb[i] = startOTest (nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount)

       for j in range (dimX * dimY):
          semb = backSemb[i][j]

          if semb > sembmax :
             sembmax  = semb   # search for maximum and position of maximum on semblance grid for given time step
             sembmaxX = latv[j]
             sembmaxY = lonv[j]

       Logfile.add ('max semblance: ' + str(sembmax) + ' at lat/lon: ' + str(sembmaxX)+','+ str (sembmaxY))
    return backSemb

# -------------------------------------------------------------------------------------------------

def execOTest (nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount) :

    f = [nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount]
    args  = Basic.floatToString (f, delim= ',')
    prog  = sys.executable + ' ' + __file__
    cmd   = prog  + ' ' + args

    Logfile.add ('--------------------------------------------', cmd)
    result = Basic.systemCmd (cmd)
    Logfile.addLines (result)
    backSemb = Basic.readVector (semb_txt)
    return backSemb

def execOTest2 () :

    for i in range (len (sys.argv)) : print sys.argv[i]

    params = Basic.stringToFloat (sys.argv[1])
    [nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount] = params
    backSemb = startOTest (int(nostat), int(nsamp), int(i), int(nstep),
                           int(dimX),   int(dimY),  mint, new_freq, int(minSampleCount), False)

    print 'backSemb = ', backSemb[0:3]
    Basic.writeVector (semb_txt, backSemb)



# -------------------------------------------------------------------------------------------------
def startOTest (nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount, isParent = True) :

    backSemb = []

    if isParent :
       backSemb = execOTest (nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount)

    else :
       trace  = Basic.readMatrix (trace_txt,  nostat, minSampleCount, '%e')
       traveltime = Basic.readMatrix (travel_txt, nostat, dimX * dimY, '%e')
       latv   = Basic.readVector (latv_txt, '%e')
       lonv   = Basic.readVector (lonv_txt, '%e')


       for j in range (dimX * dimY):
          semb  = 0
          nomin = 0
          denom = 0

          for l in range (int (nsamp)) :
             sum = 0

             for k in range (nostat) :
                relstart_samples = int ((traveltime[k][j] - mint) * new_freq + 0.5) + i * nstep

                val   =  trace[k][relstart_samples+l]
                sum   += val
                denom += (val * val)
             # endfor nostat

             nomin += sum * sum;
             semb  = nomin / (float (nostat) * denom);
          # endfor nsamp

          backSemb.append (semb)
       #endfor dimX *dimY
    #endif isParent

    return backSemb

# -------------------------------------------------------------------------------------------------

if __name__ == "__main__":
   execOTest2 ()
