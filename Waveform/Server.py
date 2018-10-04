
import sys
import os
import platform


# add local directories to import path

sys.path.append('../Common/')

import getpass
import time
import signal
from   threading import Thread
from   Queue     import Queue
import urllib2

#      Import from Common

import Basic
import Globals
import Logfile
import DataTypes

#      Constants
#
CLIENT_FLAG     = '-x'                            # options flag : if process is client


CLIENT_FLAG_SEP = '%'

#      Error code : return value of checkFkt

RETRY_IT    = -1
HAS_NO_DATA =  0
HAS_DATA    =  1

# -------------------------------------------------------------------------------------------------

class ServerCtrl(object):

    N_RETRIES   =  3 #3             # nr of retries in case of error
    N_PARALLEL  =  4 #20            # nr of parallel running clients
    WAIT_TIME   = 1.0               # Wait time in server loop in sec

    def __init__(self, nRetries = N_RETRIES, nParallel = N_PARALLEL, waitTime = WAIT_TIME,
                       printStat = False) :

        self.nRetries       = nRetries
        self.nParallel      = nParallel
        self.waitTime       = waitTime
        self.printStatistic = printStat
        self.ClientProc     = None

#endClass

# -------------------------------------------------------------------------------------------------

class ServerData(object) :

    def __init__(self) :

        self.waiting        = []
        self.running        = []
        self.finished       = []
        self.withError      = []
        self.notFound       = []
        self.withRetryFound = []
        self.hasData        = []
#endclass

# -------------------------------------------------------------------------------------------------
#
#    Server
#
# -------------------------------------------------------------------------------------------------

def _toFileName(station, ext) :


    if station.loc == '--' : loc = ''
    else :                   loc = station.loc

    return station.net + '_' + station.sta + '_' + loc + '_' + station.comp + ext

def xmlDir() :                            return os.path.join(Globals.ProtFileDir, 'xml')
def xmlFileName(station, ext) :           return os.path.join(xmlDir(), _toFileName(station, ext))

def InventoryDir() :                      return os.path.join(Globals.ProtFileDir, 'inv')
def InventoryFileName(station, ext) :     return os.path.join(InventoryDir(), _toFileName(station, ext))

def UrlDir() :                            return os.path.join(Globals.ProtFileDir, 'url')
def UrlFileName(name) :                   return os.path.join(UrlDir(), name + '.txt')

# -------------------------------------------------------------------------------------------------
#   Queue for finished client threads
#
readyQueue = Queue()

# -------------------------------------------------------------------------------------------------
#   Execute parallel system commands via threads
#
class ClientThread(Thread):

    def __init__(self, cmd, stationName, tmpFile, protFile) :

        Thread.__init__(self)
        self.cmd         = cmd          # system command
        self.execTime    = 0            # and its execution time
        self.stationName = stationName
        self.tmpFile     = tmpFile      # stdout and stderr to this file while process is running
        self.protFile    = protFile     # save it later to this file
        self.protLines   = None         # contents of the file
    # ---------------------------------------------------------------------------------------------

    def run(self):

        if os.path.isfile(self.protFile) :           # remove old protocol of client process
            os.remove(self.protFile)

        t = time.time()
        os.system(self.cmd)                          # start client synchron(thread waits here)
        self.execTime = time.time() - t               # set execution time

        if os.path.isfile(self.tmpFile) :
           os.rename(self.tmpFile, self.protFile)
           self.protLines = Basic.readTextFile(self.protFile)  # read protocol : analyse it in server loop

        readyQueue.put(self)
    # ---------------------------------------------------------------------------------------------

    def toLogfile(self, text) :
        printMsg(self.stationName, text)

# -------------------------------------------------------------------------------------------------

class ServerBase(object) :

    def __init__(self, serverName, checkFkt, ctrl=None) :

        self.serverName = serverName
        self.checkFkt   = checkFkt
        self.startFkt   = self.startClient                  # default : local fkt
        self.noDataFkt  = None                              # optional

        if ctrl == None :  self.control = ServerCtrl()
        else :             self.control = ctrl

        self.cnt        = 0
        self.startCnt   = 0
        self.d          = ServerData()

        self._init_2()
    #end

    # -------------------------------------------------------------------------------------------------

    def _init_2(self) :

        dir = Globals.ProtFileDir

        # create directory "tmp1" for this user
        #
        if not Basic.createDirectory(dir) :              Logfile.abort()

        if not Basic.createDirectory(xmlDir()) :         Logfile.abort()
        if not Basic.createDirectory(InventoryDir()) :  Logfile.abort()
        if not Basic.createDirectory(UrlDir()) :        Logfile.abort()

    # -------------------------------------------------------------------------------------------------

    def protFileName(self, station) :
        return os.path.join(Globals.ProtFileDir, self.serverName + '_' + station + '.txt')

    def tmpFileName(self, station) :
        return os.path.join(Globals.ProtFileDir, '_' + self.serverName + '_' + station + '.txt')

    def printStatistic(self) :
        intern.printStatistic_2(self.d)

     # -------------------------------------------------------------------------------------------------

    def run(self, stationList, args=None) :

        self.cnt       = 0
        self.startCnt  = 0

        if len(stationList) == 0 : return True

        Logfile.setErrorLog(True)

        intern.killRunningClients()                                 # kill running clients(of last run)
        Basic.removeFiles(Globals.ProtFileDir, self.serverName)     # remove files on dir /tmp
        Basic.removeTempFiles()

        Logfile.setErrorLog(False)

        t1 = time.time()

        try :
           self._serverLoop(stationList, args)
           endMsg = ''
           ret    = True

        except KeyboardInterrupt :
           endMsg = 'Program aborted by Control C'
           ret    = False
        #endTry

        Logfile.setErrorLog(True)
        Logfile.add(' ', endMsg)


        intern.printRunningClients()
        intern.killRunningClients()
        #endif

        Basic.removeTempFiles()
        deltaT = time.time() - t1

        if self.cnt != 0 : perStation = "%3.2f" %(deltaT / self.cnt)
        else :             perStation = '0'

        msg = 'Running ' + str(int(deltaT)) + ' sec. - ' + perStation + ' sec. per request'

        Logfile.addDelim()
        Logfile.add(' ', msg, ' ')
        Logfile.setErrorLog(False)

        return ret

    # -------------------------------------------------------------------------------------------------
    def _isClientFinished(self, stationName) :

       file = self.protFileName(stationName)

       if not os.path.isfile(file) : return False                      # client is running
       else :                         return True

    # -------------------------------------------------------------------------------------------------
    def _checkProcess(self, client, nErrors) :

        lines    = client.protLines
        execTime = client.execTime

        errorCode = self.checkFkt(client.stationName, nErrors, lines, execTime)
        return errorCode

    # -------------------------------------------------------------------------------------------------
    def startClient(self, stationName, args = None) :

        #   Special case : running client with same name(is error)
        intern.killClient(stationName)

        #   Build client args as string : stationName + ':' + arg[0] ':' arg[1] ':' ....

        if args == None : args = stationName
        else :            args = stationName + CLIENT_FLAG_SEP + args

        #  ... and command line

        file    = self.tmpFileName(stationName)
        sep     = ' '


        cmd = 'python '
        cmd += sep.join(sys.argv) + ' ' + CLIENT_FLAG + ' ' + args + ' > ' + file + ' 2>&1 '

        if self.startCnt < 3 :
           self.startCnt += 1
           Logfile.debug('cmd = ' + cmd)

        if self.control.ClientProc == None :
           protFile = self.protFileName(stationName)

           t = ClientThread(cmd, stationName, file, protFile)
           t.start()
           return t

        else :                 # ??? nur temporaer : entfernen
           oldArgs  = sys.argv
           sys.argv = cmd.split()
           self.control.ClientProc()
           sys.argv = oldArgs
        #endif
    #end

    # ---------------------------------------------------------------------------------------------

    def _serverLoop(self, stationList, args) :

        self.d.waiting = stationList
        nStations      = len(stationList)
        iMax           = min(self.control.nParallel, nStations)
        nRunning       = 0
        Logfile.add(' ')

        #
        #   Fill list of running processes with N stations
        #
        for i in range(iMax) :
            station = self.d.waiting [self.cnt]
            self.cnt += 1

            client = self.startFkt(station, args)
            client.toLogfile('Start ' + str(self.cnt) + ' / ' + str(nStations))

            nRunning += 1
            time.sleep(1.0)
        #end for

        Logfile.add(' ')

        #   Process all stations
        #
        while nRunning > 0 :
           allProcessed =(self.cnt >= len(self.d.waiting))
           t1 = time.clock()

           #
           # wait for next finished client
           #
           nWait  = 20
           client = None

           for i in range(nWait) :
               dt = 20.0

               try :
                   client   = readyQueue.get(timeout=dt)       # wait here
                   nRunning -= 1                                # next received
                   break

               except :
                   if i < 2 :         pass
                   elif i < nWait-1 : Logfile.add('Wait since ' + str((i+1) * dt) + ' sec.')
           #endfor

           if client == None :   # Nothing received since 400 sec

              if not allProcessed : Logfile.error('Continue with next task')
              else :
                  Logfile.error('No response from last ' + str(nRunning) + ' clients')
                  break             # exit main loop
           else :
              #
              #  Check finished client
              #
              station  = client.stationName
              nRetries = self.d.withError.count(station)
              ErrCode  = self._checkProcess(client, len(self.d.withError))
              timeStr  =("(%.2f sec.)  " % client.execTime)

              if ErrCode != RETRY_IT :                             # client finished ok
                 self.d.finished.append(station)

                 if ErrCode == HAS_NO_DATA :  s = ' *****'
                 else :
                    s = '(Data)'
                    self.d.hasData.append(station)
                 #endif

                 if nRetries == 0 : msg = 'Finished ' + s + timeStr
                 else :
                    msg = 'Finished ' + s + timeStr + ', retries = ' + str(nRetries)
                    self.d.withRetryFound.append(station)
                 #endif

                 client.toLogfile(msg)

                 if ErrCode == HAS_DATA : Logfile.add(' ')

              else :    # ErrCode == RETRY_IT                      # client finished with error
                 self.d.withError.append(station)

                 if self.d.withError.count(station) >= self.control.nRetries :
                    if self.control.nRetries == 1 : msg = 'finished'
                    else :                          msg = 'Retry count reached'

                    self.d.notFound.append(station)

                 else :
                    msg = 'Finished - retry later'
                    self.d.waiting.append(station)                # retry later
                 #endif

                 client.toLogfile(msg + timeStr)
              #endif ErrCode
           #endif station != None

           #  Use next station
           #
           wt = self.control.waitTime

           if wt > 0.0 :
              dt = time.clock() - t1

              if dt < wt :
                 time.sleep(wt - dt)
           #endif

           if not allProcessed :
              station = self.d.waiting [self.cnt]
              client  = self.startFkt(station, args)

              self.cnt += 1
              nRunning += 1

              startMsg =  'Start ' + str(self.cnt) + ' / ' + str(nStations)
              startMsg += '(' + str(len(self.d.hasData)) +')'

              if self.cnt > nStations :
                 startMsg += ' - retry ' + str(self.d.withError.count(station))

              client.toLogfile(startMsg)
           #endif allProcessed
        #end while

        Logfile.setErrorLog(True)
        Logfile.add(' ', 'All stations processed')
        Logfile.setErrorLog(True)

        #   print statistic

        if self.control.printStatistic :
           self.printStatistic()

        return True
#endclass  ServerBase

# -------------------------------------------------------------------------------------------------

class ClientBase(object) :

    def __init__(self, station) :
        self.station = station

    def run(self) :
        #global CLIENT_WAIT_TIME

        #try :
           # Set the signal handler(unused)
           #
           '''
           Conf   = Globals.ConfigDict
           CLIENT_WAIT_KEY  = 'duration_client'

           if CLIENT_WAIT_KEY in Conf :
              CLIENT_WAIT_TIME = int(Conf [CLIENT_WAIT_KEY])

              else :
                 signal.signal(signal.SIGALRM, intern.clientSignalhandler)
                 signal.alarm (CLIENT_WAIT_TIME)
                 pass    #hs : raus
           '''

           t1 = time.time()
           print 'call of station ', self.station
           assert self._run != None
           self._run()                               # _run() is defined in derived class

        #except :
           #print 'Client aborted by exception'
        #endtry

        #print CLIENT_TIME_MSG, int(time.time() - t1)
    #end

#endclass ClientBase

# -------------------------------------------------------------------------------------------------
'''
def clientSignalhandler(self, signum, frame):

    print 'Signal handler called with signal', signum
    print CLIENT_ABORT_MSG, CLIENT_WAIT_TIME, ' seconds.'

    time.sleep (2.0)
    sys.exit   (0)
'''
# -------------------------------------------------------------------------------------------------
# Routines for server implementation

# -------------------------------------------------------------------------------------------------
#
#   Client routine
#
CLIENT_ABORT_MSG = 'Process aborted after '
CLIENT_WAIT_TIME = 120                            # use if not set in config file
CLIENT_TIME_MSG  = 'Time = '

class ServerIntern(object):

    def __init__(self) :  dummy = 1

    # -------------------------------------------------------------------------------------------------

    def findRunningClients(self, pythonScript = sys.argv[0]):    # return filtered output of system command ps

       user = getpass.getuser()

       cmd = 'ps -AF | grep python'

       lines  = Basic.systemCmd(cmd)
       lines2 = []

       for s in lines :
           words = s.split()

           if not s.find(pythonScript) != -1 : continue
           if not s.find(CLIENT_FLAG)  != -1 : continue

           if words[0] == user : lines2.append(s)
       #endfor

       return lines2

    def printRunningClients(self) :

        stations2 = self.getStationsOfRunningClients()   # get args of all clients

        if stations2 == 0 : return

        stations  = []

        for s in stations2 :                          # extract station from client args
           sta = s.split(CLIENT_FLAG_SEP)[0]
           stations.append(sta)

        #   Print names

        names  = Basic.formatStrings(stations, '%-8s')
        sep    = ' '
        s1     = 'Running ' + str(len(names)) + ' clients :'

        if(len(names)) <= 10 : Logfile.debug(s1, sep.join(names))
        else :                   Logfile.debug(s1, sep.join(names[1:10]), sep.join(names[11:]))

    # -------------------------------------------------------------------------------------------------

    # find pid's in output of system command ps

    def _filterClientPids(self, lines) :

       pidList = []

       for s in lines :
          words = s.split()
          pid   = int(words[1])
          pidList.append(pid)
       #endfor

       return pidList

    def getClientPids(self, pythonScript = sys.argv[0]):

        lines = self.findRunningClients(pythonScript)
        return self._filterClientPids(lines)

    #  get pid of client 'name'

    def getClientPid(self, name, pythonScript = sys.argv[0]):

        lines    = self.findRunningClients(pythonScript)
        allNames = self._filterStationNames(lines)

        #print 'all names ', name, allNames

        if not name in allNames :   return -1     # not running

        pos     = allNames.index(name)
        allPids = self._filterClientPids(lines)

        return allPids [pos]

    # -------------------------------------------------------------------------------------------------
    # find station names in output of system command ps

    def _filterStationNames(self, lines) :

        stations = []

        for s in lines :
           words = s.split()
           sta   = words [len(words) - 1]
           sta   = splitClientArgs(sta)[0]

           stations.append(sta)
        #endfor

        return stations

    def getStationsOfRunningClients(self, pythonScript = sys.argv[0]):

        lines = self.findRunningClients(pythonScript)
        return  self._filterStationNames(lines)

    # -------------------------------------------------------------------------------------------------

    def killRunningClients(self, pythonScript = sys.argv[0]):                  # kill clients of last run

        isDebug = False
        pidList = self.getClientPids(pythonScript)
        nPid    = len(pidList)

        if nPid == 0 : return nPid

        Basic.killByPID(pidList)

        msg = 'kill ' + str(nPid) + ' processes'
        Logfile.debug  (msg)

    #   if Globals.isDebug() :
    #      lines = findRunningClients()
    #
    #      for s in lines : Logfile.add(s)
    #   endif

        Logfile.add(' ')
        time.sleep(3.0)

    def killClient(self, stationName) :

        pid = self.getClientPid(stationName)

        if pid == -1 : return False

        #Logfile.debug('Error : Client ' + stationName + ' running')
        pidList = []
        pidList.append(pid)
        Basic.killByPID(pidList)

        return True

    # -------------------------------------------------------------------------------------------------

    def printTable2(self, headLine, names, maxNr = -1) :

        s = headLine + ' : ' + str(len(names))

        if maxNr != -1 : s += ' / ' + str(maxNr)

        Logfile.add(' ', s, ' ')
        line    = ''

        for i in range(1, len(names)+1) :
           line +=("%-10s" % names [i-1])

           if(i % 5) == 0 : Logfile.add(line); line = ''
        #endfor

        Logfile.add(line)


    def printTable(self, headLine, names, maxNr = -1) :  # ??? noch nicht benutzt

        s = headLine + ' : ' + str(len(names))

        if maxNr != -1 : s += ' / ' + str(maxNr)

        Logfile.add(' ', s, ' ')

        sameNet = []
        line    = ''

        for i in range(0, len(names)) :
           s = names[i]
           sameNet.append(s)

           if DataTypes.isSameNetwork(sameNet[0], s) : continue

           line = ''

           for j in range(len(sameNet)) :
              line +=("%-10s" % sameNet[j])

              if j != 0 and(j % 5) == 0 :
                 Logfile.add(line)
                 line = ''
           #endfor

           sameNet = []
        #endfor

        if Logfile != '' : Logfile.add(line)

        print '--------------------------------'
        self.printTable2(headLine, names)

    # -------------------------------------------------------------------------------------------------

    def printStatistic_2(self, d) :

        finished       = list(set(DataTypes.toStationNames(d.finished)))
        withRetryFound = list(set(DataTypes.toStationNames(d.withRetryFound)))
        notFound       = list(set(DataTypes.toStationNames(d.notFound)))
        hasData        = list(set(DataTypes.toStationNames(d.hasData)))
        withError      = list(set(DataTypes.toStationNames(d.withError)))

        netFinished    = list(set(DataTypes.toNetworkNames(finished)))
        netWithData    = list(set(DataTypes.toNetworkNames(hasData)))

        Logfile.addDelim()
        Logfile.add(' ')

        nFinished    = len(finished)
        anyDataFound = False
        withoutData  = []

        for station in sorted(finished) :
            if station in hasData : anyDataFound = True
            else :                  withoutData.append(station)
        #endfor

        net = list(set(DataTypes.toNetworkNames(withoutData)))
        netWithoutData = []

        for s in net :
            if not s in netWithData : netWithoutData.append(s)

        self.printTable2('Networks with data', sorted(netWithData))

        if len(netWithoutData) == 0 :
           Logfile.add(' ','No Networks without data', ' ')
        else :
           self.printTable2('Networks without data', sorted(netWithoutData), len(netFinished))

        if len(withoutData) == 0 :
           Logfile.add('All stations with data')
        else :
           self.printTable2('Stations without data', sorted(withoutData), len(finished))

        # --------------------------------------------------------------------

        if len(withRetryFound) > 0 :
           Logfile.add(' ', 'With retry : ' + str(len(withRetryFound)), ' ')

           for station in sorted(withRetryFound) :
               if station in hasData : s = '(Data)'
               else :                  s = '       '

               printMsg(station + s, ' ', withError.count(station))
           #endfor
        #endif

        # --------------------------------------------------------------------
        Logfile.add(' ')
        return                        # ???

        if len(notFound) == 0 : Logfile.add('All stations found')
        else :
           n = str(len(notFound))
           Logfile.add('Not found : ' + n + ' after ' + str(N_RETRIES) + ' retries')

        Logfile.add(' ')

        for station in sorted(notFound) : Logfile.add(station)

        # --------------------------------------------------------------------

        Logfile.addDelim()
        return

#endclass ServerIntern

intern = ServerIntern()

# -------------------------------------------------------------------------------------------------
#
#  Global functions
#
# -------------------------------------------------------------------------------------------------

def printMsg(station, text, onlyErrorLog=False) :
    printLines(station, [text], onlyErrorLog)

def printLines(station, lines, onlyErrorLog=False) :

    lines2 = []

    for i in range(len(lines)) :
        try :
           s     = station if i == 0 else ' '
           head  =("%-15s" % s) + ' : '
           text  = lines[i]

           if text [-1] == '\n' : line = head + text[:-1]
           else :                 line = head + text

        except : line = '???'

        lines2.append(line)
    #endfor

    if onlyErrorLog : Logfile.onlyErrorLog(True)

    Logfile.addLines(lines2)

    if onlyErrorLog : Logfile.onlyErrorLog(True)

# -------------------------------------------------------------------------------------------------
def splitClientArgs(argStr) :

    sep  = CLIENT_FLAG_SEP
    args = []

    if argStr.find(sep) == -1 : args.append(argStr)
    else :                       args = argStr.split(sep)

    return args

def _joinArgs(args) :

    sep = CLIENT_FLAG_SEP
    strings = []

    for s in args : strings.append(str(s))

    return sep.join(strings)

def joinClientArgs(args, args2= None) :

    args = _joinArgs(args)
    if args2 != None : args +=(CLIENT_FLAG_SEP + _joinArgs(args2))

    return args

# -------------------------------------------------------------------------------------------------
#  Check protocol file of client
#
def checkIsTimeOut(station, traceback) :              # ??? :  Erst mal Notbehelf

    err = ''

    for i in range(len(traceback)) :
        line = traceback[i]

        #if 'ArcLinkException: Timeout waiting' in line :
           #return True

        #  getStationWaveformData.py", line 323, in make_irisrequest
        #     data = urllib.urlopen(u).read()
        #
        #  IOError: [Errno socket error] [Errno 104] Connection reset by peer
        #  IOError: [Errno socket error] [Errno 110] Connection timed out

        if 'Connection timed out' in line :
           err = 'IOError: [Errno socket error] [Errno 110] Connection timed out'; break

        if 'Connection reset by peer' in line :
           err = 'IOError: [Errno socket error] [Errno 104] Connection reset by peer'; break

        #  EOFError: telnet connection closed

        if 'telnet connection closed' in line :
           err = 'EOFError: telnet connection closed'; break

        if 'clientSignalhandler' in line :
           err = 'SignalHandler : Timeout';  break
    #endfor

    if err == '' : return False

    Logfile.setErrorLog(True)
    printMsg(station, err)
    Logfile.setErrorLog(False)

    return True
