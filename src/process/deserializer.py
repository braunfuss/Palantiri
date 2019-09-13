import sys
from palantiri.common import Basic
from collections import OrderedDict

'''
modul for deserializing pickle data from different processes
'''


def deserializeTTT(numproc, flag_rpe=False):

        L = []

        for i in range(numproc):
            if flag_rpe is True:
                data = Basic.loadDump(str(i)+'-ttt_emp.pkl')

            else:

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


def deserializeMinTMaxT(numproc, flag_rpe=False):

        L = []
        print('here')
        for i in range(numproc):
            if flag_rpe is True:
                data = Basic.loadDump('minmax-emp'+str(i)+'.pkl')
            else:
                data = Basic.loadDump('minmax-'+str(i)+'.pkl')
            if data is not None:
                L.append(data)

        mint = min([x.mint for x in L])
        maxt = max([x.maxt for x in L])

        return mint, maxt


def deserializeSembDict(numproc, flag_rpe=False):

        L = []

        for i in range(numproc):
            if flag_rpe is True:
                data = Basic.loadDump('sembDict-emp'+str(i)+'.pkl')
            else:
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


def deserializeSembMaxFile(numproc, flag_rpe=False):

    L = []

    for i in range(numproc):
        if flag_rpe is True:
            data = Basic.loadDump('sembMAX-emp'+str(i)+'.pkl')
        else:
            data = Basic.loadDump('sembMAX-'+str(i)+'.pkl')
        if data is not None:
            L.append(data)

    sembMax = OrderedDict()

    for i in L:
        for j in i.keys():
            sembMax[j] = i[j]

    return sembMax
