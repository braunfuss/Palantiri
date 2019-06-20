import os
import sys
from palantiri.common import Basic
from palantiri.common import Globals
from palantiri.common import Logfile
from collections import OrderedDict
from palantiri.process.ttt import MinTMaxT

'''
modul for deserializing pickle data from different processes
'''

def deserializeTTT(numproc):

        L = []

        for i in range(numproc):
            data = Basic.loadDump(str(i)+'-ttt.pkl')

            if data is not None:
                L.append(data)

        TTTGridMap = OrderedDict()

        for i in L:
            if sys.version_info.major >= 3:
                for j in sorted(i.keys()):
                    TTTGridMap[j] = i[j]
            else:
                for j in i.keys():
                    TTTGridMap[j] = i[j]

        return TTTGridMap

def deserializeTTT_cube(numproc):

        L = []

        for i in range(numproc):
            data = Basic.loadDump(str(i)+'-ttt.pkl')

            if data is not None:
                L.append(data)

        TTTGridMap = OrderedDict()

        for i in L:
            if sys.version_info.major >= 3:
                for j in sorted(i.keys()):
                    TTTGridMap[j] = i[j]
            else:
                for j in i.keys():
                    TTTGridMap[j] = i[j]

        return TTTGridMap

def deserializeMinTMaxT(numproc):

        L = []

        for i in range(numproc):
            data = Basic.loadDump('minmax-'+str(i)+'.pkl')

            if data is not None:
                L.append(data)

        mint = min([x.mint for x in L])
        maxt = max([x.maxt for x in L])

        return mint, maxt


def deserializeSembDict(numproc):

        L = []

        for i in range(numproc):
            data = Basic.loadDump('sembDict-'+str(i)+'.pkl')

            if data is not None:
                L.append(data)

        sembDict = OrderedDict()

        for i in L:
            if sys.version_info.major >= 3:
                for j in sorted(i.keys()):
                    sembDict[j] = i[j]
            else:
                for j in i.keys():
                    sembDict[j] = i[j]

        return sembDict


def deserializeSembMaxFile(numproc):

    L = []

    for i in range(numproc):
        data = Basic.loadDump('sembMAX-'+str(i)+'.pkl')

        if data is not None:
            L.append(data)

    sembMax = OrderedDict()

    for i in L:
        for j in i.keys():
            sembMax[j] = i[j]

    return sembMax
