import os
import platform

WINDOWS = (platform.system() == 'Windows')

from obspy.signal.trigger import recursive_sta_lta
from obspy.signal.trigger import plot_trigger as plotTrigger
from obspy.core.trace     import Trace,Stats

import numpy as np

def writeSembMaxValue(sembmaxvalue,sembmaxlat,sembmaxlon,ntimes,Config,Folder):

    fobjsembmax = open(os.path.join(Folder['semb'],'sembmaxvalue.txt'),'w')

    for i in range(ntimes):
        fobjsembmax.write ('%d %d %.20f %.2f %.2f\n'%(i,i*int(Config['step']),sembmaxvalue[i],sembmaxlat[i],sembmaxlon[i]))

    fobjsembmax.close()

def semblancestalta (sembmaxvaluevector,sembmaxlatvector,sembmaxlonvector):

    data=np.array(sembmaxvaluevector, dtype=np.float64)

    #stats = Stats()                  #hs
    #stats.network= 'BW'          #hs
    #stats['station'] = 'MANZ'        #hs

    tr = Trace(data,header=None)

    sta = 0.5
    lta = 4
    cft = recursive_sta_lta (tr, int(sta * tr.stats.sampling_rate), int(lta * tr.stats.sampling_rate))

    #print cft

    thrOn  = 0.5
    thrOff = 1.5
    plotTrigger(tr, cft, thrOn, thrOff)
