
import os
import sys
import platform

WINDOWS = (platform.system() == 'Windows')

# add local directories to import path  

sys.path.append ('../tools/') 
sys.path.append ('../Common/')                              
           
from optparse import OptionParser 

import fnmatch
import fileinput
import re
import logging
import multiprocessing
import time
import subprocess

import getpass                          
from   os  import stat

import datetime
import signal

import urllib
from   ConfigParser import SafeConfigParser   

import obspy.core
from   obspy.arclink.client import Client 
from   obspy.core.util      import NamedTemporaryFile
from   obspy.core           import Stream, read
from   obspy.mseed.core     import isMSEED

if not WINDOWS :
   from   pwd  import getpwuid   

import config

# -------------------------------------------------------------------------------------------------

#      Import from Common

import Globals                     # Own global data
import Basic                       # Own module with basic functions
import Logfile                     # Implements logfile
import ConfigFile                  # Semantic of config file entries
import Debug
import DataTypes                   # Common data types

#      Import from Waveform

from   Version  import  VERSION_STRING

import KeyFile
import DataDir
import Server

# -------------------------------------------------------------------------------------------------
          
logger = logging.getLogger('event_station_diff')
logger.setLevel (logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)

logger.addHandler(ch)

options = None
args    = None
_mail   = None
_key_search_folder = None
_mseed_search_folder = None
_channel = None
_loc = None

# ------------------------------------------------------------------------------------------------
# <a href="http://www.hersteller.de" target="_blank">hier</a>

def saveUrl (station, url) :

    return # ???

    if Globals.isDebug :  
       file = Server.UrlFileName ('data_url')

       if url == None : Basic.writeTextFile (file, list (' \n'))
       else :   
          lines = []
          lines.append (station + ' : ' + url + '\n')

          Basic.appendToFile (file, lines)
    #endif

# ------------------------------------------------------------------------------------------------

def make_time (begin_time,duration):
   
    b_time_1     = obspy.core.utcdatetime.UTCDateTime (begin_time,iso8601=True) 
    geofon_begin = b_time_1.formatArcLink()

    e_time      = b_time_1 + float(duration) * 60
    geofon_end  = e_time.formatArcLink()
    
    b_time_iris = obspy.core.utcdatetime.UTCDateTime (begin_time,iso8601=True) 
    iris_begin  = b_time_iris.formatIRISWebService()

    e_time_iris = b_time_iris + float(duration) * 60
    iris_end    = e_time_iris.formatIRISWebService()

    dict_time = {'g_begin':geofon_begin, 'g_end':geofon_end,
                 'i_begin':iris_begin,   'i_end':iris_end,
                 'obspy_begin':b_time_1, 'obspy_end':e_time}

    return dict_time

# -------------------------------------------------------------------------------------------------

def keyfolder2List (p):

    Logfile.red ('Parsing KEYFILE Structure')
    
    stationfilespath=p
    L = []

    pattern = 'station_' + "[A-Za-z]"
    content = os.listdir(stationfilespath)

    for i in content:
        if re.match(pattern, i):
            name = i[8:]
            L.append(name)

    L_sauber = list(set(L))
    L_sauber.sort()

    K = []

    for i in L_sauber:
        if options.stationversion == 'acq':
            f    = 'station_'+i
            path = os.path.join (stationfilespath, f)

            for line in fileinput.input(path) :
                if "PACKAGES" in line:
                    if fnmatch.fnmatch(line, '*acquisition*') :
                        K.append(i)

        if options.stationversion == 'all':  K.append(i)
        if options.stationversion == '':     Logfile.add ('set station method')
        
    Logfile.red ('FOUND: --> ' + str(len(K)) + ' KEYFILES')

    return K

# -------------------------------------------------------------------------------------------------

def sdsList(p):

    L = []
    Logfile.red ('Parsing SDS Structure')
    
    for root,dirs,files in os.walk(p):
        for i in files:
            line = str.split(i,'.')
            val= line[0]+'_'+line[1]
            L.append(val)

    L_sauber = list(set(L))

    Logfile.red ('FOUND: --> ' + str (len (L_sauber)) + ' STATIONS')   
    return L_sauber

# -------------------------------------------------------------------------------------------------

def cmpSdsKey(K,S):

    Logfile.red ('Parsing for Missing Station')

    t = list (set(K).difference(set(S)))
    t = set  (t)
    Logfile.add (str (len(t)) + ' STATIONS are missing')

    return t

# -------------------------------------------------------------------------------------------------

def to_sds(output):
    return                                                           # hs ???

    cmd = 'scart -I '+output+' '+_mseed_search_folder
    p   = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()
    logger.info('\033[31munpacking %s to SDS\n\n \033[0m' % (cmd))

# -------------------------------------------------------------------------------------------------

def create_dir (dir):
  if os.access (dir, os.W_OK) :  return True

  try :    os.makedirs(dir)
  except : return False

  return True

# -------------------------------------------------------------------------------------------------

def multiplex (filename):

    st = read (filename)

    for trace in st:
        path = os.path.join (_mseed_search_folder, str (trace.stats.starttime.year),
                             trace.stats.network,  trace.stats.station, ('%s.D') % (trace.stats.channel))
        create_dir (path)

        diff      = trace.stats.endtime.day - trace.stats.starttime.day
        rangediff = diff+1

        if diff == 0:
            d         = obspy.core.utcdatetime.UTCDateTime (trace.stats.starttime)
            finalname = DataDir.filename (trace, d)
            filepath  = os.path.join (path,finalname)
            trace.write (filepath,format='MSEED',reclen=512)

        elif diff >= 1:
            for i in xrange (rangediff):
                mult = 60 * 60 * 24 * i

                if i == 0:
                    s1 = trace.stats.starttime
                    e1 = obspy.core.utcdatetime.UTCDateTime (year=s1.year, month=s1.month, day=s1.day,hour = 23, 
                                                             minute=59,second=59, microsecond=999999)
                else:
                    s1 = trace.stats.starttime+mult
                    s1 = obspy.core.utcdatetime.UTCDateTime (year=s1.year, month=s1.month, day=s1.day,hour = 0, 
                                                             minute=0, microsecond=00)
                    e1 = obspy.core.utcdatetime.UTCDateTime (year=s1.year, month=s1.month, day=s1.day,hour = 23,
                                                             minute=59,second=59, microsecond=999999)
                if i == diff:
                    e1 = trace.stats.endtime

                d         = obspy.core.utcdatetime.UTCDateTime (s1)
                finalname = DataDir.filename (trace, d)
                filepath  = os.path.join (path,finalname)

                tr        = trace.slice  (s1, e1)
                tr.write (filepath,format='MSEED',reclen=512)
                print s1,' to ',e1,jd,finalname,filepath
        #endif        
    #endfor

    #hs+
    #  "remove" unused channels
    #
    channels  = []                                                          

    for trace in st : channels.append (trace.stats.channel)      # get all channels from trace list

    channnels = sorted (channels)
    group     = [['BHE','BHN','BHZ'], ['BH1','BH2','BHZ'], ['HHE','HHN','HHZ']]  # use only these
    used      = []
    
    for chn in channels :                                                     # build channel sets
        if len (used) > 0 : break

        for i in range(3) :
            if chn in group[i] :  
               used = group[i]        # use this set 
               break;
    #endfor

    trace = st[0]
    path  = os.path.join (_mseed_search_folder, str (trace.stats.starttime.year),
                          trace.stats.network,  trace.stats.station)

    for chn in channels :
        if chn in used : continue

        # replace unused channels (files) with magic number '4711'

        dir = os.path.join (path, ('%s.D') % (chn))

        if os.path.exists (dir) :
           files = os.listdir (dir)
           line  = '4711'

           for file in files :
               Basic.writeTextFile (os.path.join (dir,file), line)
               #Logfile.add (file) 
    #endfor
    #hs-
       
# -------------------------------------------------------------------------------------------------

def proof_file_v3 (filename,component):
    
    if len(filename) == 0 : return 0          #hs

    print 'PROOF FILE V3 ',filename

    size       = 0
    streamList = []
    t          = filename[0]
    savename   = '{n}.{s}.{l}.ms'.format (n=t.stats.network, s=t.stats.station,l=t.stats.location)

    filename.write (savename, format='MSEED', reclen=512)
          
    try:
       st = read (savename)

       for i in st : streamList.append (i.stats.channel)
 
       print 'CHANNEL: ',len(streamList),len(component)

       if len (set(streamList)) == len(component) :
          #to_sds    (savename)                      #hs
          #multiplex (filename)                      #hs
          multiplex (savename)                       #hs
          size = os.path.getsize (savename)

       else:                                          
          #size = 0                                  #hs
          multiplex (savename)                       #hs
          size = os.path.getsize (savename)          #hs

    except :
       size = 0
       Logfile.exception ('proof_file_v3')           

    os.remove (savename)
    return size

# -------------------------------------------------------------------------------------------------

def make_arcObspyrequest (station,begin,end,channel,loc,pwdkeys):
    
    print 'INOBSPY ',pwdkeys,type(pwdkeys)
    print 'Obspy ARCLINKREQUEST FOR STATION '

    size = 0
    rstation = station.split('_')

    if loc == '--': loc = ''

#   client = Client (user=_mail,dcid_keys=pwdkeys,timeout = 5)        #hs - Lutz
    client = Client (user=_mail,dcid_keys=pwdkeys,timeout = 120)      #hs

    st     = Stream ()
    
    #for i in channel:
    #    try :    st += client.getWaveform (rstation[0], rstation[1], loc, i, begin, end)
    #    except : continue         # ArclinkException : No data available
    #endfor

    try :    st += client.getWaveform (rstation[0], rstation[1], loc, '*', begin, end)  #hs v2
    except : return 0

    print 'INOBS ',rstation,channel    
    streamList = []

    for i in st : 
        streamList.append (i.stats.channel)
    
    #if len (set (streamList)) == len (channel) : size = proof_file_v3 (st, channel) #hs
    #else:                                        size = 0                           #hs
    
    size = proof_file_v3 (st, channel)            #hs
    return size

# -------------------------------------------------------------------------------------------------------------------------------

def make_irisrequest (station,begin,end,channel,loc):
    
        print 'IRISREQUEST FOR STATION ' , station

        net,sta = station.split('_')
        url = 'http://service.iris.edu/fdsnws/dataselect/1/query?'
        st  = Stream()

        for i in channel:

            parameter = urllib.urlencode({
                 'net': net,
                 'sta': sta,
                 'loc': loc,
                #'cha': i,                               #hs v2
                 'starttime': begin,
                 'endtime': end,
                 'nodata': '404',
            })

            u     = ('%s%s')%(url,parameter)
            saveUrl (station, u)
        
            #data = urllib.urlopen(u).read()                          #hs
            try :    data = urllib.urlopen(u).read()                  #hs
            except : continue

            #tf  = NamedTemporaryFile()                               #hs
            tf   = NamedTemporaryFile (suffix = 'xtmp')               #hs
            tf.write(data)
            tf.seek(0)

            t = isMSEED(tf.name)

            if t == True:
                st += read(tf.name, 'MSEED')

            tf.close()
            os.remove(tf.name)
            break                                                #hs v2 

#       if t == True:     size = proof_file_v3 (st,channel)      #hs
        if len(st) > 0 :  size = proof_file_v3 (st,channel)      #hs
        else:             size = 0

        print 'SIZE: ----> ',size  
        return size

# -------------------------------------------------------------------------------------------------------------------------------

def write_statistic (year,day):

    K = keyfolder2List(_key_search_folder)
    S = sdsList(_mseed_search_folder)
    MISS = cmpSdsKey(K,S)

    fobjm = open (os.path.join(os.getcwd(),'missing-station-'+str(year)+'-'+str(day)+'.dat'),'w')
    fobjm.write  ('\n'.join(MISS))
    fobjm.close  ()
    
    fobja = open (os.path.join(os.getcwd(),'available-station-'+str(year)+'-'+str(day)+'.dat'),'w')
    fobja.write  ('\n'.join(S))
    fobja.close  ()


# -------------------------------------------------------------------------------------------------------------------------------

def make_meta(evpath):
    
    #d = obspy.core.utcdatetime.UTCDateTime(stime)
    #jd = "%03d" % d.julday
    #cmd='python ev_meta_mt4.py -p '+_mseed_search_folder+' -y '+str(d.year)+' -d '+str(jd)
    #write_statistic(d.year,jd)
    #print cmd
    
    evpath = evpath.split('/')[-1]
    logger.info('\033[31mNEXT PROCESSING STEP: \n\n 1) python arraytool.py getmeta {evdirectory} \n\n\033[0m'.format(evdirectory=evpath))
    
    #p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #p.wait()
    #os.system(cmd)


# -------------------------------------------------------------------------------------------------------------------------------

def stationRequest_old (station,timeDict,counter,pwdkeys,maxstation):
    
        for i in _loc:
            for j in _channel:
                print 'loc, channel = ', i,j
               
                size = make_irisrequest(station,timeDict['i_begin'],timeDict['i_end'],j,i)

                if size != 0: print  'DATA FROM IRIS FOR STATION'; return
                else:
                    size = make_arcObspyrequest (station,timeDict['obspy_begin'],
                                                         timeDict['obspy_end'],j,i,pwdkeys)

                    if size != 0: print 'DATA FROM WEBDC FOR STATION'; return
                    else:         print 'NO DATA FROM WEBDC OR IRIS FOR STATION AND FINISH'
                    #endif
                #endif

                print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n'
            #endfor   
        #endfor       


NO_DATA_FLAG = 'NO DATA FROM WEBDC OR IRIS FOR STATION'

def stationRequest (station, timeDict, counter, pwdkeys, maxstation):
    
        isIris = KeyFile.isIRIS (options.keyfolder, fullName = station)
        
        for i in _loc:
           j = '*'
           print 'loc, channel = ', i,j

           try :
              if isIris :
                 if make_irisrequest (station, timeDict['i_begin'],timeDict['i_end'],j,i) != 0 :
                    print  'DATA FROM IRIS FOR STATION'
                    return
              else :
                 if make_arcObspyrequest (station, timeDict['obspy_begin'] ,timeDict['obspy_end'], j, i, pwdkeys) != 0 :
                    print 'DATA FROM WEBDC FOR STATION'
                    return

           except : return Logfile.exception ('stationRequest')
        #endfor       
 
        print NO_DATA_FLAG

# -------------------------------------------------------------------------------------------------

def globalConf():

    cDict  = {}
    parser = SafeConfigParser()
    parser.read (os.path.join ('..', 'global.conf'))

    for section_name in parser.sections():
        for name, value in parser.items(section_name) : cDict[name] = value

    return cDict

# -------------------------------------------------------------------------------------------------

def checkNoData (station, traceback) :              # ??? :  Erst mal Notbehelf

    ret = False
    err = None

    for i in range (len (traceback)) :
        line = traceback[i]

        # in make_irisrequest
        #  if len(st) > 0 :  size = proof_file_v3 (st,channel)      #hs
        #    in proof_file_v3
        #     filename.write (savename, format='MSEED', reclen=512)
        #     in write
        #       writeFormat(self, filename, **kwargs)
        #          in writeMSEED
        #            (1.0 / trace.stats.sampling_rate * HPTMODULUS) % 100 != 0:
        #               ZeroDivisionError: float division by zero

        if 'trace.stats.sampling_rate' in line :
           err = 'Sampling rate is zero'; 
           ret = True
           break

        if 'HTTPError' in line :
           err = 'HTTP Error 404: Not Found'
           ret = True
           break
    #endfor

    if err != None :
       Logfile.setErrorLog (True)
       Server.printMsg (station, 'Error : ' + err)
       Logfile.setErrorLog (False)
    #endif

    return ret

# -------------------------------------------------------------------------------------------------

def checkProcessError (station, nErrors, lines, execTime) :   # ??? execTime

    errCode = Server.HAS_NO_DATA
    errCode = Server.RETRY_IT

    keyfile  = KeyFile.KeyFileObj (options.keyfolder, fullName = station)
    sta      = keyfile.read ()

    if sta == None : site = None
    else :           site = sta.site + ' (' + sta.provider + ')'

    for lineNr in range (len (lines)) :
       line  = lines [lineNr]
       isEnd = False
       s     = ' '

       # UserWarning: MAX_REQUESTS exceeded - breaking current request loop

       if 'MAX_REQUESTS' in line : 
          errCode = Server.RETRY_IT
          s =  'UserWarning: MAX_REQUESTS exceeded - breaking current request loop'
          s += ' (' + str(nErrors) + ')'
          #isEnd = True
 
       elif 'deprecated' in line : s = ' '             # ignore ObsPyDeprecation Warning   #15.7.2016

       elif Logfile.MSG_TOKEN         in line :  s = line
       elif 'python'                  in line :  s = line

       elif Server.CLIENT_ABORT_MSG   in line :
          errCode = Server.RETRY_IT
          s       = line

       elif line.startswith ('DATA FROM') :
          s     = line
          isEnd = True

       # PROOF FILE V3  3 Trace(s) in Stream:
       # AK.FID..BHE | 2014-04-03T02:43:14.008400Z - 2014-04-03T04:43:13.988400Z | 50.0 Hz, 360000 samples
       # AK.FID..BHN | 2014-04-03T02:43:14.008400Z - 2014-04-03T04:43:13.988400Z | 50.0 Hz, 360000 samples
       # AK.FID..BHZ | 2014-04-03T02:43:14.008400Z - 2014-04-03T04:43:13.988400Z | 50.0 Hz, 360000 samples
       # CHANNEL:  3 3
       # SIZE: ---->  961024

       elif 'PROOF FILE' in line  :                 # station has data
          errCode = Server.HAS_DATA
          sn      = []

          if site != None : sn.append (site)

          for i in range (0,100) :
              if lineNr+i >= len (lines) : break

              if 'CHANNEL'        in lines [lineNr+i] : break
              if '#### Exception' in lines [lineNr+i] : break

              sn.append (lines [lineNr+i])
          #endfor

          lineNr += i
          Server.printMsg   (' ', ' ')
          Server.printLines (station, sn)
       #endif

       elif 'Traceback' in line :                        # client aborted with traceback info
          sn = []

          for i in range (0,300) :
              if lineNr+i >= len (lines) : break         #10.12.2015
              if 'KeyboardInterrupt' in lines [lineNr+i] : sn = []; break
             #if lineNr+i >= len (lines) : break         #10.12.2015

              sn.append (lines [lineNr+i])
          #endfor

          if Server.checkIsTimeOut (station, sn) :       # Traceback shows timeout
             Logfile.error ('Retry access later')
             errCode = Server.RETRY_IT

          elif checkNoData (station, sn) :               # Traceback shows 'no data'
             errCode = Server.HAS_NO_DATA
            #errCode = Server.RETRY_IT
             return errCode

          else :                                         # Traceback --> log
             Server.printLines (station, sn, onlyErrorLog=True)

          isEnd = True
       #endif

       if s != ' ' : Server.printMsg (station, s)
       if isEnd    : break
    #endwhile

    return errCode

# -------------------------------------------------------------------------------------------------

def initWaitingList (parameter) :

    K       = keyfolder2List (parameter.keyfolder) 
    S       = sdsList        (parameter.sdsfolder)
    MISS    = cmpSdsKey (K,S)
    MISS    = list (set (MISS))
    waiting = []

    n = len (MISS)

    for i in range (n) : waiting.append (MISS[i])
    return waiting

# --------------------------------------------------------------------------------------------------

def init (isClient) :

    Globals.isClient = isClient
    Debug.init ()

    if not isClient :
       if not Logfile.init (startMsg = VERSION_STRING) : return False

    return Globals.init () 

# --------------------------------------------------------------------------------------------------

def checkConfigFile (conf) :

    mail          = ConfigFile.mail
    pwd           = ConfigFile.pwd
    keyfilefolder = ConfigFile.keyfilefolder
    duration      = ConfigFile.duration

    keyList = [mail, pwd, keyfilefolder, duration]  
    return ConfigFile.checkKeys (conf, keyList)

# --------------------------------------------------------------------------------------------------
    
SERVER_NAME = 'data'

def startIrisServer (stations) :

    ctrl = Server.ServerCtrl (nRetries = 1, nParallel=10, waitTime=0.1, printStat=False)
    srv  = Server.ServerBase (SERVER_NAME, checkProcessError, ctrl)
    #if WINDOWS : srv.control.ClientProc = MainProc

    return srv.run (stations)


def startGeofonServer (stations) :
          
    ctrl = Server.ServerCtrl (nRetries = 2, nParallel=10, waitTime=1.0, printStat=False)
    srv  = Server.ServerBase (SERVER_NAME, checkProcessError, ctrl)
    srv.control = ctrl
    #if WINDOWS : srv.control.ClientProc = MainProc

    return srv.run (stations)

# -------------------------------------------------------------------------------------------------
#
#   Client routine 
#
class WaveformClient (Server.ClientBase) :
    
    def __init__ (self, options, pwdkeys) :
        Server.ClientBase.__init__ (self, options.station)

        self.options = options
        self.pwdkeys = pwdkeys

    def _run (self) :  # called direct from ClientBase

        time_d = make_time (self.options.time, self.options.duration)
        stationRequest (self.options.station ,time_d, 0, self.pwdkeys, 0)

# --------------------------------------------------------------------------------------------------

def run_parallel (options, pwdDict) :

    if options.station :                                       # Client part
       if not init (True) : return False

       clt = WaveformClient (options, pwdDict)
       clt.run ()
       return True

    else :                                                    # Server part
       if not init (False) : return False

       keyfileDir = os.path.join (Globals.EventDir(), options.keyfolder)

       if not Basic.checkExistsDir (keyfileDir) : 
          return False                                 # Cannot find directory

       #   Build station list
       #
       stationList = sorted (initWaitingList (options))

       if len (stationList) == 0 :
          return Logfile.error ('No stations found')

       #  Check program version
       #
       if not KeyFile.checkVersion (keyfileDir, fullName = stationList[0]) :
          return False

       #  Run server(s)
       #
       saveUrl (' ', None)        # init debug service
       network = options.network

       mask       = KeyFile.getIrisMask (None, stations=stationList) 
       irisList   = Basic.selectStrings (stationList, mask)
       geofonList = Basic.selectStrings (stationList, Basic.Not (mask))

       if not network or network == 'iris' :
          if not startIrisServer (irisList) :
             return True                             # aborted with ctrl c
       #endif
  
       if not network or network == 'geofon' :
          if not startGeofonServer (geofonList) :
             return True                             # aborted with ctrl c
      #endif  

       if network and network != 'iris' and network != 'geofon' :
          if not KeyFile.isNetwork (network) :
             return Logfile.error ('Illegal network name <' + network + '>', 
                                   'Network not found in directory ' + Globals.KeyfileFolder())

          list2 = DataTypes.selectNetwork (irisList, network)     # search in iris list

          if len(list2) > 0 :
             startIrisServer (list2)
             return True                         
         
          list2 = DataTypes.selectNetwork (geofonList, network)   # search in geofon list

          if len(list2) > 0 :
             startGeofonServer (list2)
             return True                         
       #endif

#hs- ----------------------------------------------------------------------------------------------

def main (args):
    
    parser = OptionParser (usage="\npython %prog -t 2009-12-31T12:23:34 -d 5 -m SDS -k key -s all/acq")

    parser.add_option ("-t","--time",      type="string", dest="time",      help="time")
    parser.add_option ("-d","--duration",  type="string", dest="duration",  help="duration in min")
    parser.add_option ("-m","--sdsfolder", type="string", dest="sdsfolder", help="sdsfolder")
    parser.add_option ("-k","--keyfolder", type="string", dest="keyfolder", help="keyfolder")
    parser.add_option ("-s","--station",   type="string", dest="stationversion", help="stationversion")
    parser.add_option ("-f","--evpath",    type="string", dest="eventpath", help="eventpath")
    parser.add_option ("-x","--dummy",     type="string", dest="station",   help="dummy")    #hs : client flag
    parser.add_option ("-n","--dummy2",    type="string", dest="network",   help="dummy2")   #hs : single network

    return parser.parse_args(args)

# -------------------------------------------------------------------------------------------------

def MainProc () :

    global options,args
    global _mail,_key_search_folder,_mseed_search_folder,_channel,_loc

    options,args = main (sys.argv)

    Basic.checkExistsDir (options.eventpath, isAbort=True)
    Globals.setEventDir  (options.eventpath)

    C      = config.Config (options.eventpath)
    Origin = C.parseConfig ('origin')
    Conf   = globalConf()

    checkConfigFile (Conf)

    _mail = Conf['mail']
    
    pwdDict = {}

    for i in Conf['pwd'].split(','):
        pwdDict[i.split(':')[0]]=i.split(':')[1]

    #print pwdDict
    
    options.time           = Origin ['time']
    options.duration       = int (Conf['duration'])
    options.sdsfolder      = os.path.join (options.eventpath,'data')
    options.keyfolder      = os.path.join (options.eventpath, Conf['keyfilefolder'])
    options.stationversion = 'all'
    
    _key_search_folder   = options.keyfolder
    _mseed_search_folder = options.sdsfolder
    
    channel1 = ['BHE','BHN','BHZ']
    channel2 = ['BH1','BH2','BHZ']
    channel3 = ['HHE','HHN','HHZ']

    _channel = [channel1,channel2,channel3]
    _loc     = ['--','10','00','11']

    run_parallel (options, pwdDict)
                                         
    if not options.station :                          # server
       Logfile.showLabel ('Programm finished')

# -------------------------------------------------------------------------------------------------
#  Aufruf : python  tools/getStationWaveformData.py  -f  <event dir>
#
if __name__ == "__main__" :

    MainProc()
