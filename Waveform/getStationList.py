
import os
import sys
import platform

WINDOWS = (platform.system() == 'Windows')

# add local directories to import path

sys.path.append ('../tools/') 
sys.path.append ('../Common/')                              
 
from   optparse import OptionParser
import logging
import time
import multiprocessing
import fnmatch

import xml.dom.minidom
import csv

from datetime import datetime

import  config
import  urllib
import  urllib2   
from    lxml     import etree   

from ConfigParser import SafeConfigParser

from obspy.arclink.client     import Client   
from obspy.core.util          import locations2degrees
from obspy.core.utcdatetime   import UTCDateTime

import cPickle as pickle 

if WINDOWS :
   import xlrd              # for convertion xls to cvs
#endif

#       Import from Common

import  Basic
import  Globals
import  Logfile
import  Debug
import  Server
import  ConfigFile                  # Semantic of config file entries
from    ConfigFile import ConfigObj

#        Import from Waveform

from    Version  import  VERSION_STRING

import  KeyFile
import  WebDC                      # Arclink client

HAS_DATA = 'has data'

def printMsg (station, text=' ') :
    try    :  Server.printMsg (Logfile.MSG_TOKEN + station, text)
    except :  dummy = 1

# -------------------------------------------------------------------------------------------------

class Station(object):
    '''
    Class for storing station information
    
    :type net: str
    :param net: Network of station
    :type station: str
    :param station: Name of station
    :type lat: str
    :param lat: Latitude of Station
    :type lon: str
    :param lon: Longitude of Station
    :type elev: str
    :param elev: Elevation of Station
    :type site: str
    :param site: Location of Station
    :type atime: str
    :param atime: Archive time of IRIS data stream , optional
    :type rtime: str
    :param rtime: Realtime data stream of IRIS, optional
    :type start: str
    :param start: Archive time of WEBDC data stream,optional
    :type end: str
    :param end: Realtime data stream of WEBDC, optional
    :type provider: str
    :param provider: Data provider of the station, optional
    
    '''
    
    def __init__(self, net, station, lat, lon, elev, site,
                 atime=None,rtime=None,start=None,end=None,provider=None):
        '''
        Initializes the station object
        '''
        self.net   = net;    self.station = station; self.lat      = lat;     self.lon = lon
        self.elev  = elev;   self.atime   = atime;   self.rtime    = rtime;   self.site = site
        self.start = start;  self.end     = end;     self.provider = provider

    def __str__(self):
        return ('%s.%s')%(self.net,self.station)

# -------------------------------------------------------------------------------------------------
#
#
class NetworkList(object):
    '''
    Class of network list client to retrieve network list of WEBDC and IRIS
    
    :type otime: str
    :param otime: Event time 
    :type elat: str
    :param elat: Event latitude
    :type elon: str
    :param elon: Event longitude
    :type minDist: str
    :param minDist: Minimum distance from station to event
    :type maxDist:  str
    :param maxDist: Maximum distance from station to event
    :type duration: str
    :param duration: duration from event time on where data should be available for the station 
    :type mail: str
    :param mail: Mail adress of user to get network information from WEBDC , optional 
    :type blacklist: list
    :param blacklist: List of Networks who will be blacklisted in the search
    '''
    
    def __init__(self,otime,elat,elon,minDist,maxDist,duration,
                 mail='ehlert@geo.uni-potsdam.de',blacklist=[]):
        '''
        Initializes the NetworkList object
        '''
        self.otime   = otime;   self.elat      = elat;     self.elon      = elon
        self.minDist = minDist; self.maxDist   = maxDist;  self.duration  = duration
        self.mail    = mail;    self.blacklist = blacklist

# -------------------------------------------------------------------------------------------------

    def getIRISList(self):
        '''
        Returns list of available IRIS networks
        '''
        inetworks = self._listIrisNetworks()
        inetworks = list (set (inetworks))
        
        self._filterStationsBlacklist (inetworks)
        
        return sorted (inetworks)

# -------------------------------------------------------------------------------------------------

    def getWEBDCList(self):
        '''
        Returns list of available WEBDC networks
        '''

        gnetworks  = self._listGeofonNetworks()   
        gnetworks2 = self._listGeofonNetworks2()

        if WINDOWS : gnetworks3 = []  # ??? noch fuer Windows einbauen
        else :
           t = UTCDateTime(self.otime)        
           gnetworks3 = self._getGeofonNetworks (t-10, t+10) 
        
        gnetworks.extend (gnetworks2)  
        gnetworks.extend (gnetworks3) 
        gnetworks = list (set(gnetworks))
        
        self._filterStationsBlacklist (gnetworks)

        return sorted (gnetworks)
# -------------------------------------------------------------------------------------------------
    
    def _listIrisNetworks(self):
        #Downloads available IRIS networks, removes duplicates and returns them as list
        
        data = []

       #URL  = 'http://www.iris.edu/dms/nodes/dmc/services/network-codes/?type=csv' #9.12.2015
        URL  = 'http://www.fdsn.org/networks/?type=csv'                             #9.12.2015
        s    = 'download latest IRIS permanent network tables : '

        if Basic.existsHTML_Page (URL, s, withComment = True) :
           Logfile.add (' ', s, URL)
           datareader = csv.reader (urllib2.urlopen (URL))

           for row in datareader : data.append (row[0])
        #endif

        URL = 'http://www.iris.edu/SeismiQuery/bin/tempNetsExcel.php'   
        s   = 'download latest IRIS temporary network tables : '

        if Basic.existsHTML_Page (URL, withComment = True) :
           Logfile.add (' ', s, URL)

           tempnetname = Globals.TempFileName ('allTempNetstemp.xls')
           tmpcsv      = Globals.TempFileName ('allTempNetstemp.csv')

           u = urllib2.urlopen (URL)

           localFile = open (tempnetname, 'w')
           localFile.write  (u.read())
           localFile.close  ()

           try :
              if WINDOWS : 
                 wb = xlrd.open_workbook (tempnetname)
                 
              else :       
                 os.system (('in2csv  %s > %s') % (tempnetname,tmpcsv))
                 datareader = csv.reader (open (tmpcsv, 'rb'), delimiter=",")

              for row in datareader:
                  if len (row[0]) == 2 : data.append (row[0])

           except :
              Logfile.error ('Cannot convert to cvs file')

           if os.path.isfile (tmpcsv) :       os.remove (tmpcsv)
           if os.path.isfile (tempnetname) :  os.remove (tempnetname)
        #endif

        data = sorted (list (set(data)))
        return data

# -------------------------------------------------------------------------------------------------
    
    def _listGeofonNetworks2 (self):
        '''
        Download available WEBDC networks from EIDA
        '''
    
        L   = []
        URL = 'http://www.orfeus-eu.org/Data-info/eida-station-overview.txt'
        s   = 'download latest EIDA network tables :'

        if not Basic.existsHTML_Page (URL, withComment = True) : 
           return L

        Logfile.add (' ', s, URL)

        #nettxt = urllib2.urlopen (URL)                                              #hs
        tmpFile = Globals.TempFileName ('eida_index.txt')                            #hs
        nettxt  = Basic.readURL (URL, tmpFile)                                       #hs

        for line in nettxt:
            line = line.split()

            if line[5] == 'OPEN': L.append (line[0])
        #endfor

        L = sorted (list (set(L)))
        return L
# -------------------------------------------------------------------------------------------------

    def _listGeofonNetworks_old (self):                     #hs : Code von Lutz - nicht mehr benutzt
        '''
        Download available networks via Geofon kml file
        '''

        URL = 'http://geofon.gfz-potsdam.de/waveform/archive/kml.php'
        Logfile.add (' ','download latest GEOFON network tables :')

        L = []
        kml_file = urllib2.urlopen (URL)
        baum     = etree.parse (kml_file)
        tag_dict = baum.getroot()

        for eintrag in tag_dict.getchildren():
            for t in eintrag.getchildren():
                for i in t:
                    if i.tag == '{http://www.opengis.net/kml/2.2}name' :  L.append(i.text)
        return L

    def _listGeofonNetworks (self):                         #hs : new routine : replaces ..._old
        '''
        Download available networks via Geofon kml file 
        '''      
        return WebDC.listNetworks ()

# -------------------------------------------------------------------------------------------------
   
    def _getGeofonNetworks (self,start,end):
        '''
        Return dictionary of available networks via Arclink
        
        :type start: obspy.core.utcdatetime.UTCDateTime
        :param start: Start date and time
        :type end: obspy.core.utcdatetime.UTCDateTime
        :param end: End date and time
        '''
        return WebDC.getNetworks (self.mail, start, end)

# -------------------------------------------------------------------------------------------------  
    
    def _filterStationsBlacklist (self,NList):
        '''
        Delete blacklisted networks from network list
        
        :type NList: list
        :param Nlist: List of Networks to delete blacklisted Networks
        '''

        for i in self.blacklist:
            try:     del NList[NList.index(i)]
            except:  continue

# -------------------------------------------------------------------------------------------------

RETRY_ACCESS = 'Retry access'

def getNetworkInventory (network):
        '''
        Retrieve all stations from one WEBDC network
        
        :type network: str
        :param network: name of network to search for station
        '''
        
        inv = WebDC.getNetworkInventory (network, 'ehlert@geo.uni-potsdam.de')

        if inv == None : print RETRY_ACCESS

        return inv

# -------------------------------------------------------------------------------------------------

def parseInventory (Dict):
        '''
        Parses Network dictionary from WEBDC networks to retrieve available stations
        
        :type Dict: dictionary
        :param Dict: network dictionary with all station information
        '''

        StationList = []
        prov = KeyFile.PROV_WEB_DC                   #hs

        for i in Dict.iterkeys():
            t = i.split('.')

            if len(t) == 2:
                net    = t[0]
                sta    = t[1]
                lat    = Dict[i]['latitude']
                lon    = Dict[i]['longitude'] 
                elev   = Dict[i]['elevation']
                site   = Dict[i]['description']
                tstart = Dict[i]['start']
                tend   = Dict[i]['end']

#               newSta = Station (net,sta,lat,lon,elev,site,start=tstart,end=tend)                #hs
                newSta = Station (net,sta,lat,lon,elev,site,start=tstart,end=tend, provider=prov) #hs
                StationList.append (newSta)
        #endfor

        return StationList

def parseInventory_new (Dict):
        '''
        Parses Network dictionary from WEBDC networks to retrieve available stations
        
        :type Dict: dictionary
        :param Dict: network dictionary with all station information
        '''
        return WebDC.parseInventory (Dict)

# -------------------------------------------------------------------------------------------------

def buildGeofonMsg (stationobject, text=None) :

    start = str (stationobject.start).split ('T')[0]
    end   = str (stationobject.end).split   ('T')[0]
    now   = str (UTCDateTime (datetime.now())).split ('T')[0]

    if end.strip() == now.strip() : end = '           '                           

    s = start + ' - ' + end + '  ' + toStr (stationobject.site)
    
    if text != None : s += (' (' + text +')')

    return s

# -------------------------------------------------------------------------------------------------
    
def filterStationsTimeGeofon (stationList, parameter):
        '''
        Filter stations via time attribute
        
        :type stationList: list
        :param stationList: list of stationobjects
        :type parameter: list
        :param parameter: list of parameter used filter routines
        '''

        t0      = UTCDateTime (parameter['time'])
        msgList = []
        errList = []
        G       = []

        for stationobject in stationList:
            if stationobject.end == None:
                stationobject.end = UTCDateTime (datetime.now())

            if UTCDateTime (stationobject.end) > t0 and UTCDateTime (stationobject.start) < t0 :
               msgList.append (buildGeofonMsg (stationobject))
               G.append (stationobject)

            else : msgList.append (buildGeofonMsg (stationobject))
        #endfor

        #   Print stations with time interval

        if len (msgList) > 0 :
           printMsg (' ')
           printMsg ('GEOFON network ')

           for i in range (len (stationList)) :
               printMsg (str (stationList[i]), msgList[i])

           print HAS_DATA
        #endfor

        return G
# -------------------------------------------------------------------------------------------------   
    
def create_dir (directory):
    '''
    Creates recursively directorys with access control
    
    :type directory: str
    :param directory: name of the path to be created
    '''
    if os.access (directory, os.W_OK): return True

    try:
        os.makedirs(directory)
        return True

    except: return False

# -------------------------------------------------------------------------------------------------

def toStr (s) :
    try : return str (s)
    except : return '???'

# -------------------------------------------------------------------------------------------------

def filterStationsDistance (stationList, parameter):
        '''
        Filters stationslist via distance attribute
        
        :type stationList: list
        :param stationList: list of stationobjects
        :type parameter: list
        :param parameter: list of parameter used filter routines
        '''
    
        D = []

        for stationobject in stationList:
            sdelta = locations2degrees (parameter['elat'], parameter['elon'], 
                                        float (stationobject.lat), float (stationobject.lon))

            if sdelta > parameter['minDist'] and sdelta < parameter['maxDist']:
               D.append(stationobject)
        #endfor

        return D
# -------------------------------------------------------------------------------------------------

def buildIRISMsg (stationobject, start, end, text) :

    start = ('%s-%s-%s') % (start[0],start[1],start[2])
    end   = ('%s-%s-%s') % (end[0],  end[1],  end[2])
    now   = str (UTCDateTime (datetime.now())).split ('T')[0]

    if end.strip() == now.strip() : end = '           '                           

    s = start + ' - ' + end + '  ' + toStr (stationobject.site)
    
    if text != None : s += (' (' + text +')')

    return s


def filterStationsTimeIRIS (stationList, parameter):
    '''
    Filters stations via time attribute
    
    :type stationList: list
    :param stationList: list of stationobjects
    :type parameter: list
    :param parameter: list of parameter used filter routines
    '''

    t0      = UTCDateTime (parameter['time'])
    msgList = []
    F       = []

    for stationobject in stationList:
        t = stationobject.atime.split ('-')
        z = stationobject.rtime.split ('-')

        if len(t) == 2:
            start = t[0].split('/')
            end   = t[1].split('/')
            s     = ('%s-%s-%sT00:00:00.0Z') % (start[0],start[1],start[2])
            e     = ('%s-%s-%sT00:00:00.0Z') % (end[0],  end[1],  end[2])

            if  UTCDateTime(e) > t0 and UTCDateTime(s) < t0 :
                msgList.append (buildIRISMsg (stationobject, start, end, 'ARDATA'))
                F.append (stationobject)
            else :
                msgList.append (buildIRISMsg (stationobject, start, end, 'ARDATA - outside'))

        if len(z) == 2:
            start = z[0].split('/')
            end   = z[1].split('/')
            s     = ('%s-%s-%sT00:00:00.0Z') % (start[0],start[1],start[2])
            e     = ('%s-%s-%sT00:00:00.0Z') % (end[0],  end[1],  end[2])

            if  UTCDateTime(e) >= t0 and UTCDateTime(s) <= t0 :
                msgList.append (buildIRISMsg (stationobject, start, end, 'RTDATA'))
                F.append (stationobject)
            else :
                msgList.append (buildIRISMsg (stationobject, start, end, 'RTDATA - outside'))
    #endfor

    if len (msgList) > 0 :
       printMsg (' ')
       printMsg ('IRIS network ')

       for i in range (len (stationList)) :
           try    : printMsg (str (stationList[i]), msgList [i])
           except : continue

       print HAS_DATA                                                 #hs ???
    #endfor

    return F

# -------------------------------------------------------------------------------------------------

def parseXML (request):
    '''
    Parse XML representation of IRIS networks to retrieve station information
    
    :type request: str
    :param request: name of iris network
    '''
    StationList = []

    baseurl = 'http://www.iris.edu/cgi-bin/xmlstationinfo/'
    request = baseurl+request
    prov    = KeyFile.PROV_IRIS

    try:
            doc = xml.dom.minidom.parse (urllib2.urlopen(request))

            for node in doc.getElementsByTagName("station"):
                net   = str (node.getAttribute("net"))
                sta   = node.getAttribute("sta")
                lat   = node.getAttribute("lat")
                lon   = node.getAttribute("lon") 
                elev  = node.getAttribute("elev")
                atime = node.getAttribute("ardata")
                rtime = node.getAttribute("rtdata")
                site  = node.getAttribute("site")

#               newSta = Station (net,sta,lat,lon,elev,site,atime,rtime)                  #hs
                newSta = Station (net,sta,lat,lon,elev,site,atime,rtime, provider=prov)   #hs
                StationList.append (newSta)
    except:
        Logfile.error ("NOT PARSABLE " + request)

    return StationList
# -------------------------------------------------------------------------------------------------

def geofonMt (geofonnetwork,pid,parameter):
        '''
        Function to retrieve network information for WEBDC networks
        
        :type geofonnetwork: list
        :param geofonnetwork: list of networks from webdc
        :type pid: int
        :param pid: process id for pickle file name
        :type parameter: list
        :param parameter: list of parameter used filter routines
        '''
        
        SDict = getNetworkInventory (geofonnetwork)
        T     = []

        if SDict != None :   
           S     = parseInventory (SDict)
           D     = filterStationsDistance   (S,parameter)
           T     = filterStationsTimeGeofon (D,parameter)
        
        return T

# -------------------------------------------------------------------------------------------------

def irisMt (irisnetwork,pid,parameter):
    '''
    Function to retrieve network information for IRIS networks
    
    :type irisnetwork: list
    :param irisnetwork: list of networks from iris
    :type pid: int
    :param pid: process id for pickle file name
    :type parameter: list
    :param parameter: list of parameter used filter routines
    '''
    
    S = parseXML (irisnetwork)
    D = filterStationsDistance (S,parameter)
    T = filterStationsTimeIRIS (D,parameter)

    return T
  
# -------------------------------------------------------------------------------------------------

def checkProcessError (station, nErrors, lines, execTime = None) :

    errCode = Server.RETRY_IT
    #errCode = Server.HAS_NO_DATA

    #  Print messages from client
    #
    t1 = Logfile.MSG_TOKEN
    n  = 2 * len(t1)

    for s in lines :
        if t1 in s : Logfile.add (s[n:-1])

    Logfile.add (' ')

    #   Check if errors in log
    #
    for lineNr in range (len (lines)) :
       line  = lines [lineNr]
       isEnd = False
       s     = ' '

       # UserWarning: MAX_REQUESTS exceeded - breaking current request loop -> retry access

       if 'MAX_REQUESTS' in line : 
          errCode = Server.RETRY_IT
          s =  'UserWarning: MAX_REQUESTS exceeded - breaking current request loop'
          s += ' (' + str(nErrors) + ')'
 
       elif 'deprecated' in line : s = ' '             # ignore ObsPyDeprecation Warning   #15.7.2016

       elif 'python' in line : s = line

       elif 'Signal handler' in line :
          errCode = Server.RETRY_IT; s = line

       elif Server.CLIENT_ABORT_MSG in line :
          errCode = Server.RETRY_IT;  s = line

       elif HAS_DATA in line  :                   # station has data
          errCode = Server.HAS_DATA
          isEnd   = True

       elif 'Traceback' in line :                 # client aborted with traceback info
          sn = []

          for i in range (0,300) :
              if lineNr+i >= len (lines) : break         #10.12.2015
              if 'KeyboardInterrupt' in lines [lineNr+i] : sn = []; break
             #if lineNr+i >= len (lines) : break         #10.12.2015

              sn.append (lines [lineNr+i])
          #endfor

          if Server.checkIsTimeOut (station, sn) :       # Traceback shows timeout
         #if True :
             Logfile.error ('Retry access later')
             errCode = Server.RETRY_IT

          else :                                         # Traceback --> log
             Server.printLines (station, sn, onlyErrorLog=True)

          isEnd = True
       #endif

       if s != ' ' : Server.printMsg (station, s)
       if isEnd    : break
    #endwhile

    return errCode

# --------------------------------------------------------------------------------------------------
 
def init (options) :

    isClient = (options.args != None)

    Globals.isClient = isClient
    Debug.init ()

    if not isClient :
       if not Logfile.init (startMsg = VERSION_STRING) : return False

       Basic.checkExistsDir (options.evpath, isAbort=True)

    Globals.setEventDir  (options.evpath)

    return Globals.init ()

# --------------------------------------------------------------------------------------------------

def checkConfigFile (conf) :

    mail          = ConfigFile.mail
    pwd           = ConfigFile.pwd
    duration      = ConfigFile.duration
    minDist       = ConfigFile.mindist
    maxDist       = ConfigFile.maxdist
    blackList     = ConfigFile.blacklist

    keyList = [mail, pwd, duration, minDist, maxDist]  
    ConfigFile.checkKeys (conf, keyList)

    keyList = [blackList]
    ConfigFile.checkKeys (conf, keyList, optional=True)

# --------------------------------------------------------------------------------------------------

def run_parallel (options) :

    '''
    Starts station search procedure
    
    :type options: instance
    :param options: parameter to initialize the networklist class
    '''
    isClient = (options.args != None)

    if not init (options) : 
       return False
 
    if isClient :                                       #        Run client

       clt = StationListClient (options)
       clt.run ()
       return True 

    else :                                              #         Run server
       #    Create directory for clients
       #
       clientDir = os.path.join (options.evpath, 'keyfiles-' + str (time.time()))

       Logfile.add ('Create keyfile directory ', clientDir, ' ')  
       create_dir  (clientDir)

       #  Build network list
       #
       C      = config.Config (options.evpath)
       Origin = C.parseConfig ('origin')
       Conf   = Globals.ConfigDict
    
       checkConfigFile (Conf)  
      
       globalCfg = ConfigObj (dict = Conf)
       originCfg = ConfigObj (dict = Origin)
    
       ot       = originCfg.Time()                          # str   (Origin['time'])
       elat     = originCfg.lat()                           # Origin['lat']
       elon     = originCfg.lon()                           # Origin['lon']

       minDist  = globalCfg.Distance ('mindist')            # Conf  ['mindist']
       maxDist  = globalCfg.Distance ('maxdist')            # Conf  ['maxdist']
       duration = globalCfg.Duration ()                     # Conf  ['duration']

       paramList = [ot, maxDist,minDist,elat,elon]

       BL = []

       if 'blacklist' in Conf :
          K  = (Conf ['blacklist']).split(',')
          BL = ['# Network Code']
          BL.extend(K)
       
       T = NetworkList (ot,elat,elon,minDist,maxDist,duration, blacklist = BL,mail=Conf['mail'])

       SERVER_NAME = 'network'

       #         Handle Iris networks
       #
       inetworks = T.getIRISList()   
      #inetworks = ['BF']  
    
       if len (inetworks) == 0 :  Logfile.error ('No iris networks found')
       else : 
          args = Server.joinClientArgs ([IRIS_TAG, clientDir], paramList)
          ctrl = Server.ServerCtrl (nRetries = 1, nParallel=1, waitTime=1.0, printStat=False)    
          srv  = Server.ServerBase (SERVER_NAME, checkProcessError, ctrl)

          #if WINDOWS : srv.control.ClientProc = MainProc

          if not srv.run (inetworks, args) : return False
       #endif
       
       #       Handle Geofon networks
       #
       gnetworks = T.getWEBDCList() 
      #gnetworks = ['AF']
      #gnetworks = ['FR']
      #gnetworks = []
   
       if len (gnetworks) == 0 :
          Logfile.error ('No geofon networks found')

       else :
          #     Access network infos now from Geofo
          #
          args = Server.joinClientArgs ([GEOFON_TAG, clientDir], paramList)
          ctrl = Server.ServerCtrl (nRetries = 4, nParallel=1, waitTime=2.0, printStat=False)    
          srv  = Server.ServerBase (SERVER_NAME, checkProcessError, ctrl)

          #if WINDOWS : srv.control.ClientProc = MainProc

          if not srv.run (gnetworks, args) : return False
       #endif
    
       #    Print statistic
       
       nIres  = len (inetworks)
       nWebDC = len (gnetworks)
       nAll   = nIres + nWebDC

       if nIres  != 0 : Logfile.add (' ', 'Processed ' + str(nIres)  + ' IRES networks')
       if nWebDC != 0 : Logfile.add (     'Processed ' + str(nWebDC) + ' WEBDC networks')

       if nAll == 0 : return Logfile.error ('No networks found')

       if   nIres  == 0 : err = 'No IRIS network found'
       elif nWebDC == 0 : err = 'No WEBDC network found'
       else :             err = None

       if err != None : Logfile.add (err)
       
       # showNextStep
       #
       evpath        = options.evpath.split('/')[-1]
       keyfoldername = clientDir.split('/')[-1]

       Logfile.add (' ', 'NEXT PROCESSING STEP:', ' ')
       Logfile.add ('   1) change keyfolder value in global.conf to ' + keyfoldername)
       Logfile.add ('   2) python arraytool.py getdata ' + evpath, ' ')

       return True
    #endif  # end of server

# -------------------------------------------------------------------------------------------------
#
#   Client routine 
#
IRIS_TAG   = 'i'
GEOFON_TAG = 'g_'

class StationListClient (Server.ClientBase) :
    
    def __init__ (self, options) :

        print 'Start client ', options.args
        args  = Server.splitClientArgs (options.args)

        self.net     = args[0]
        self.tag     = args[1]
        self.dirname = args[2]
        ot           = args[3]
        maxDist      = float (args[4])
        minDist      = float (args[5])
        elat         = float (args[6])
        elon         = float (args[7])

        self.param = {'time':ot,   'maxDist':maxDist, 'minDist':minDist,
                      'elat':elat, 'elon':elon }

        Server.ClientBase.__init__ (self, self.net)

    def _run (self) :  # called direct from ClientBase

        if   self.tag == IRIS_TAG :   stationList = irisMt   (self.net, -1, self.param)
        elif self.tag == GEOFON_TAG : stationList = geofonMt (self.net, -1, self.param)
        else :                        assert False

        stationList = list (set (stationList))
   
        for stat in stationList :
            # Write dummy Key Files with reduced station information
            #
            keyfile = KeyFile.KeyFileObj (self.dirname, stat.net, stat.station)
            keyfile.write (stat)
            #keyfile.read  ()
        #endfor

        if len(stationList) > 0 :  print HAS_DATA 

#endclass StationListClient

# -------------------------------------------------------------------------------------------------

def main (args):
    '''
    Parse commandline arguments
    
    :type args: list
    :param args: List of commandline arguments to parse for the parameter of the program
    '''
    parser = OptionParser (usage="%prog -f eventpath ")

    parser.add_option ("-f", "--evpath", type="string", dest="evpath", help="evpath")
    parser.add_option ("-x", "--dummy",  type="string", dest="args",   help="dummy")   #hs : client flag

    return parser.parse_args(args)

# --------------------------------------------------------------------------------------------------

def MainProc () :

     options,args = main (sys.argv)
     run_parallel (options)

     if not options.args :                          # server
        Logfile.showLabel ('Program finished')

if __name__ == '__main__':

   MainProc ()
