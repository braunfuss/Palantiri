import os
import sys
import time

from threading import Thread


# add local directories to import path
sys.path.append ('../Common/')

import numpy as np          #    Import from NumPy

import Logfile              #    Import from common
import Basic


trace_txt  = 'trace.txt'
travel_txt = 'travel.txt'
latv_txt   = 'latv.txt'
lonv_txt   = 'lonv.txt'
semb_txt   = 'semb.txt'


def startC_Code (nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_freq, minSampleCount) :

        f    = [nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_freq, minSampleCount]
        args = Basic.floatToString (f, delim= ',')
        path =  os.path.dirname (__file__)

        prog = os.path.join (path,'OTest')

        os.system (prog  + ' ' + args)

# -------------------------------------------------------------------------------------------------

class MyThread (Thread):

    def __init__ (self, nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount):
        Thread.__init__(self)

        self.nostat    = nostat
        self.nsamp     = nsamp
        self.i         = i
        self.nstep     = nstep
        self.dimX      = dimX
        self.dimY      = dimY
        self.mint      = mint
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
               latv_1, lonv_1, traveltime_1, trace_1) :

   #if USE_C_CODE  :
    USE_C_CODE = True
    if USE_C_CODE == True :
       return otestSeriell (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_frequence,
                            minSampleCount, latv_1, lonv_1, traveltime_1, trace_1)
    else :
       return otest_py     (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_frequence,
                            minSampleCount, latv_1, lonv_1, traveltime_1, trace_1)

def otest_py  (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_frequence, minSampleCount,
               latv_1, lonv_1, traveltime_1, trace_1) :

    trace      = toMatrix (trace_1, minSampleCount)
    traveltime = toMatrix (traveltime_1, dimX * dimY)
    latv       = latv_1.tolist()
    lonv       = lonv_1.tolist()

    '''
    Basic.writeMatrix (trace_txt,  trace, nostat, minSampleCount, '%e')
    Basic.writeMatrix (travel_txt, traveltime, nostat, dimX * dimY, '%e')
    Basic.writeVector (latv_txt,   latv, '%e')
    Basic.writeVector (lonv_txt,   lonv, '%e')
    '''

    backSemb = np.ndarray (shape=(ntimes, dimX*dimY), dtype=float)
    print "lob"
    for i in range (ntimes) :
        #  loop over grid points
        sembmax = 0; sembmaxX = 0; sembmaxY = 0

        for j in range (dimX * dimY):
            semb = 0; nomin = 0; denom = 0

            for l in range (int (nsamp)) :
               sum = 0
               for k in range (nostat) :
                   relstart = int ((traveltime[k][j] - mint) * new_frequence + 0.5) + i * nstep
                   sum     += trace[k][relstart + l]
                   denom   += trace[k][relstart + l] * trace[k][relstart + l]

            nomin += sum * sum

            x    = latv[j]
            y    = lonv[j]
            semb = nomin / (float (nostat) * denom)

            #print 'sembout ', x,y,semb
            backSemb[i][j] = semb

            if semb > sembmax :
               sembmax  = semb   # search for maximum and position of maximum on semblance
                                 # grid for given time step
               sembmaxX = latv[j]
               sembmaxY = lonv[j]


        Logfile.add ('max semblance: ' + str(sembmax) + ' at lat/lon: ' +
                     str(sembmaxX)+','+ str (sembmaxY))

    return backSemb

# -------------------------------------------------------------------------------------------------

def otestSeriell (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_freq, minSampleCount,
                  latv_1, lonv_1, traveltime_1, trace_1) :

    trace      = toMatrix (trace_1, minSampleCount)
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
    print "noppe3"
    return backSemb

# -------------------------------------------------------------------------------------------------

def otestPar (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_freq, minSampleCount,
              latv_1, lonv_1, traveltime_1, trace_1) :

    trace      = toMatrix (trace_1, minSampleCount)
    traveltime = toMatrix (traveltime_1, dimX * dimY)
    latv       = latv_1.tolist()
    lonv       = lonv_1.tolist()

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
       sembmax     = 0
       sembmaxX    = 0
       sembmaxY    = 0
       backSemb[i] = startOTest (nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount)

       for j in range (dimX * dimY):
          semb = backSemb[i][j]

          if semb > sembmax :
             sembmax  = semb   # search for maximum and position of maximum on semblance grid for given time step
             sembmaxX = latv[j]
             sembmaxY = lonv[j]
          #endif
       #endfor dimX *dimY

       Logfile.add ('max semblance: ' + str(sembmax) + ' at lat/lon: ' + str(sembmaxX)+','+ str (sembmaxY))
    #endfor ntimes
    print "nope4"
    return backSemb

# -------------------------------------------------------------------------------------------------

def execOTest (nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount) :

    f     = [nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount]
    args  = Basic.floatToString (f, delim= ',')
    prog  = sys.executable + ' ' + __file__
    cmd   = prog  + ' ' + args

    Logfile.add ('--------------------------------------------', cmd)
    result = Basic.systemCmd (cmd)
    Logfile.addLines (result)
    print "nope2"
    backSemb = Basic.readVector (semb_txt)
    return backSemb

def execOTest2 () :

    for i in range (len (sys.argv)) : print sys.argv[i]

    params = Basic.stringToFloat (sys.argv[1])
    [nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount] = params

    print 'Enter startOTest'

    backSemb = startOTest (int(nostat), int(nsamp), int(i), int(nstep),
                           int(dimX),   int(dimY),  mint, new_freq, int(minSampleCount), False)

    print 'backSemb = ', backSemb[0:3]
    Basic.writeVector (semb_txt, backSemb)

    print 'Leave startOTest'


# -------------------------------------------------------------------------------------------------
def startOTest (nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount, isParent = True) :

    backSemb = []

    if isParent :
       backSemb = execOTest (nostat, nsamp, i, nstep, dimX,dimY, mint, new_freq, minSampleCount)

    else :
       trace      = Basic.readMatrix (trace_txt,  nostat, minSampleCount, '%e')
       traveltime = Basic.readMatrix (travel_txt, nostat, dimX * dimY, '%e')
       latv       = Basic.readVector (latv_txt, '%e')
       lonv       = Basic.readVector (lonv_txt, '%e')

       print
       #  loop over grid points
       #
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
   print "yes"
   execOTest2 ()
