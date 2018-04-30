
import os
import platform

WINDOWS = (platform.system() == 'Windows')

from   obspy.arclink.client import Client

#       Import from Common

import  Basic
import  Globals
import  Logfile

# -------------------------------------------------------------------------------------------------
# used by getStationList.py :

def _getFromCatalog (network) :

    return None   

    # spaeter ???

    dir      = os.path.join (Globals.EventDir(), "keyfiles_catalog")
    files    = os.listdir (dir)
    selected = []

    for file in files :
        if file.startwith ('station_' + network + '_') :
           selected.append (file)
    #endfor



def getNetworkInventory (network, user1):
        '''
        Retrieve all stations from one WEBDC network
        
        :type network: str
        :param network: name of network to search for station
        '''
        
        inv  = None
        isOk = False

        for i in range(5) :
            try :
               client = Client (user=user1)                                            #hs
               #client  = Client (user=user1, timeout=100)                             #hs 
               inv     = client.getInventory (network, '*', '*', '*', restricted=None,permanent=None)
               isOk    = True
               break

            except : 
             #time.sleep (1.0); continue
              return None
        #endfor

        if not isOk :  return None   # return _getFromCatalog (network)
        else :         return inv


# used by ev_meta.py :           

def getInventory (station, usermail_1, key_1) :

    client = Client (user=usermail_1, dcid_keys=key_1)
    inv    = client.getInventory (station.net, station.sta, loc, station.comp, instruments=True)

    return inv

# -------------------------------------------------------------------------------------------------
# used by getStationList.py

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

# -------------------------------------------------------------------------------------------------
# used by getStationList.py

def getNetworks (user1, start, end) :
    '''
    Return dictionary of available networks via Arclink
        
    :type start: obspy.core.utcdatetime.UTCDateTime
    :param start: Start date and time
    :type end: obspy.core.utcdatetime.UTCDateTime
    :param end: End date and time
    '''

    for i in range (3) :
       Logfile.add (' ','Waiting for dictionary of available geofon networks (via Arclink) ....')

       try :
          L        = []
          client  = Client (user = user1)                           #hs
          #client = Client (user = user1, timeout=20)               #hs
          t       = client.getNetworks (start,end)
    
          for index,i in enumerate (t.iterkeys()):
             z = i.split('.')
             L.append (z[0])

          L = list (set(L))
          break

       except :
          Logfile.exception ('getGeofonNetworks')
          Logfile.exception ('Retry access')
          L = []
          continue
    #endfor

    return L

# -------------------------------------------------------------------------------------------------
# used by getStationList.py

def listNetworks_old ():                     #hs : Code von Lutz - nicht mehr benutzt
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

# -------------------------------------------------------------------------------------------------

def listNetworks ():                         #hs : new routine : replaces ..._old
    '''
    Download available networks via Geofon kml file
    '''

    L    = []
    URL  = 'http://geofon.gfz-potsdam.de/waveform/archive/index.php?type=p'     # permanent
    #URL = 'http://geofon.gfz-potsdam.de/waveform/archive/index.php?type=all'   # all
    s    = 'download latest GEOFON network tables :'

    if not Basic.existsHTML_Page (URL, s, withComment = True) : 
       Logfile.error ('Cannot find url :', URL)
       return L

    Logfile.add (' ',s, URL)

    tmpFile = os.path.join  (Globals.ProtFileDir, 'geofon_index.txt') 
    lines   = Basic.readURL (URL, tmpFile)

    for line in lines :
        word = 'network='

        if 'station.php' in line and word in line :              
           pos     = line.find (word) + len(word)
           network = line [pos: pos+2]
           #Logfile.add (network)
           L.append (network)
        #endfor

    L = list (set(L))
    return L

