
'''
modul for deserializing pickle data from different processes
'''

import os
import sys

sys.path.append('../tools/')
sys.path.append('../Common/')

if sys.version_info.major >= 3:
    import _pickle as cPickle
else:
    import cPickle as pickle
#       Import from common

import  Basic
import  Globals
import  Logfile
from collections import OrderedDict

#       Import from Process

from ttt import MinTMaxT

# -------------------------------------------------------------------------------------------------

def deserializeTTT(numproc):

        L = []

        for i in range(numproc):
            data = Basic.loadDump(str(i)+'-ttt.pkl')

            if data != None: L.append(data)
        #endfor

        TTTGridMap = OrderedDict()

        for i in L:
            if sys.version_info.major >= 3:
                for j in sorted(i.keys()) : TTTGridMap[j] = i[j]
            else:
                for j in i.iterkeys(): TTTGridMap[j] = i[j]

        return TTTGridMap

# -------------------------------------------------------------------------------------------------

def deserializeMinTMaxT(numproc):

        L = []

        for i in range(numproc):
            data = Basic.loadDump('minmax-'+str(i)+'.pkl')

            if data != None: L.append(data)
        #endfor

        mint = min([x.mint for x in L])
        maxt = max([x.maxt for x in L])

        return mint,maxt

# -------------------------------------------------------------------------------------------------

def deserializeSembDict(numproc):

        L = []

        for i in range(numproc):
            data = Basic.loadDump('sembDict-'+str(i)+'.pkl')

            if data != None: L.append(data)
        #endfor

        sembDict = OrderedDict()

        for i in L:
            if sys.version_info.major >= 3:
                for j in sorted(i.keys()) : sembDict[j] = i[j]
            else:
                for j in i.iterkeys(): sembDict[j] = i[j]

        return sembDict


def deserializeSembMaxFile(numproc):

    L = []

    for i in range(numproc):
        data = Basic.loadDump('sembMAX-'+str(i)+'.pkl')

        if data != None:  L.append(data)
    #endfor

    sembMax = OrderedDict()

    for i in L:
        for j in i.iterkeys(): sembMax[j] = i[j]

    return sembMax
