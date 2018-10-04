
import os
import sys


sys.path.append('../Common/')

import  Basic
import  Globals
import  Logfile
import  Debug

import  DataTypes


PROV_WEB_DC = 'WEB_DC'
PROV_IRIS   = 'IRIS'
PROV_UNDEF  = '???'

class KeyFileObj(object):


   def __init__(self, dirName=None, net=None, station=None, fullName=None) :

       if dirName == None : dir = Globals.KeyfileFolder()
       else :               dir = dirName

       if fullName != None :
          net     = DataTypes.toNetwork(fullName)
          station = DataTypes.toStation(fullName)

       self.dirName  = dir
       self.fullName = net + '_' + station
       self.key      = None

   # ----------------------------------------------------------------------------------------------

   def dirName (self) :        return self.dirName
   def fullName(self) :        return self.fullName

   # ----------------------------------------------------------------------------------------------

   def _keyfileName(self, net, station) :

       fname = 'station_' + str(net) + '_' + str(station)
       return os.path.join(self.dirName, fname)

   def _error(self, text) :
       Logfile.error(self.fullName + ',' + self.key + ' : ' + text)

   def _floatError(self, text) :
       _error(self, text)
       return 0.0

   def _Float(self, s, minVal=None, maxVal=None) :

       try :   val = float(s[1:-1])
       except :
           _error(self, 'Cannot convert to float ' + s)
           return 0.0

       if minVal != None :
          if val < minVal :
             return _floatError(self, str(val) + ' - Range error(< ' + str(minVal))

       if maxVal != None :
          if val > maxVal :
             return _floatError(self, str(val) + ' - Range error(> ' + str(minVal))

       return val


   def _String(self, s) :
       try :    return str(s)
       except :
          _error(self, 'Illegal string')
          return '???'

   # ----------------------------------------------------------------------------------------------

   def read(self):

       net      = DataTypes.toNetwork(self.fullName)
       station  = DataTypes.toStation(self.fullName)
       fname    = self._keyfileName  (net, station)

       if not os.path.isfile(fname) : return None

       lines    = Basic.readTextFile (fname)

       if len(lines) == 0 : return None

       sta = DataTypes.Station(net, station, loc='???',comp='???')

       try :
          END_FLAG = 'PACKAGES'
          endFound = False

          for i in range(len(lines)) :
             line = lines [i].strip()
             #print 'line= ', line

             w      = line.split('=')
             key    = self._String(w[0])
             _g_Key = key
             val    = w[1]

             if   key ==  'KEY_VERSION'     : dummy        = self._Float(val)                #  0 KEY_VERSION='2.5'
             elif key ==  'STAT_DESC'       : sta.site     = val                              #  1 STAT_DESC='Ganja, Azerbaijan'
             elif key ==  'LATITUDE'        : sta.lat      = self._Float (val, -90.0, 90.0)  #  2 LATITUDE='40.6519'
             elif key ==  'LONGITUDE'       : sta.lon      = self._Float (val,-180.0,360.0)  #  3 LONGITUDE='46.3297'
             elif key ==  'ELEVATION'       : sta.ele      = self._Float (val)               #  4 ELEVATION='560.0'
             elif key ==  'DATALOGGER'      : dummy        = self._String(val)               #  5 DATALOGGER='Q380-M'
             elif key ==  'DATALOGGER_SN'   : dummy        = self._String(val)               #  6 DATALOGGER_SN='xxxx'
             elif key ==  'SEISMOMETER1'    : dummy        = self._String(val)               #  7 SEISMOMETER1='STS-2N'
             elif key ==  'SEISMOMETER_SN1' : dummy        = self._String(val)               #  8 SEISMOMETER_SN1='yyyy'
             elif key ==  'GAIN_MULT1'      : dummy        = self._Float (val)               #  9 GAIN_MULT1='1.0'
             elif key ==  'SAMPLING1'       : dummy        = self._String(val)               # 10 SAMPLING1='20/40/80/100'
             elif key ==  'DEPTH1'          : dummy        = self._Float (val)               # 11 DEPTH1='0.0'
             elif key ==  'SEISMOMETER2'    : dummy        = self._String(val)               # 12 SEISMOMETER2=''
             elif key ==  'SEISMOMETER_SN2' : dummy        = self._String(val)               # 13 SEISMOMETER_SN2=''
             elif key ==  'GAIN_MULT2'      : dummy        = self._String(val)               # 14 GAIN_MULT2=''
             elif key ==  'SAMPLING2'       : dummy        = self._String(val)               # 15 SAMPLING2=''
             elif key ==  'DEPTH2'          : dummy        = self._String(val)               # 16 DEPTH2=''
             elif key ==  'START_DATE'      : dummy        = self._String(val)               # 17 START_DATE='1980/001'
             elif key ==  'CONFIGURED'      : dummy        = self._String(val)               # 18 CONFIGURED='yes'
             elif key ==  'PACKAGES'        : sta.provider = self._String(val) [1:-1]        # 19 PACKAGES='WEB_DC'

             else                           : #self._error('Invalid key ' + key)
                Logfile.error('Invalid key ' + key)

             if key == END_FLAG :
                endFound = True
                break
          #endfor

       except :
          Logfile.exception('readKeyFile', fname)

       if not endFound : Logfile.error(self.fullName+ ' : keyfile cut')

       return sta

   # -------------------------------------------------------------------------------------------------

   def toStr1(self,s) :                  #hs : str() in writeKeyFile geht manchmal nicht
        try : return str(s)
        except : return '???'

   def write(self, stationobject):
        '''
        Write dummy Key Files with reduced station information
        '''
        fname = self._keyfileName(stationobject.net, stationobject.station)
        fobj  = open(fname, 'w')

        try:
           msg ='''KEY_VERSION='2.5'
STAT_DESC=\'''' + self.toStr1(stationobject.site) + '''\'
LATITUDE=\''''  + self.toStr1(stationobject.lat)  + '''\'
LONGITUDE=\'''' + self.toStr1(stationobject.lon)  + '''\'
ELEVATION=\'''' + self.toStr1(stationobject.elev) + '''\'
DATALOGGER='Q380-M'
DATALOGGER_SN='xxxx'
SEISMOMETER1='STS-2N'
SEISMOMETER_SN1='yyyy'
GAIN_MULT1='1.0'
SAMPLING1='20/40/80/100'
DEPTH1='0.0'
SEISMOMETER2=''
SEISMOMETER_SN2=''
GAIN_MULT2=''
SAMPLING2=''
DEPTH2=''
START_DATE='1980/001'
CONFIGURED='yes'
PACKAGES=\'''' + self.toStr1(stationobject.provider) + '''\'
       '''
           fobj.write(msg)

        except :
           Logfile.exception('writeKeyFile', fname)

        fobj.close()


#end class KeyfileObj
# --------------------------------------------------------------------------------------------------

def getNetworks(dirName=None) :

    if dirName == None : dir = Globals.KeyfileFolder()
    else :               dir = dirName

    files    = os.listdir(dir)
    networks = []

    for s in files :
        word = str.split(s,'_')

        if len(word) == 3 and word[0] == 'station' :
           networks.append(word[1])
    #endfor

    return sorted(set(networks))

def isNetwork(network, dirName = None) :

    assert network != None
    return network in getNetworks(dirName)

# --------------------------------------------------------------------------------------------------

#def getStation(net=None, station=None, fullName=None) :  ???

    #file = KeyFileObj(dirName=None, net,station, fullName)
    #sta  = file.read()
    #return sta

def getProvider(dirName=None, net=None, station=None, fullName=None) :

    file = KeyFileObj(dirName, net,station, fullName)
    sta  = file.read()

    if sta == None : return None
    else :           return sta.provider


def getSite(stationName) :

    keyfolder = Globals.KeyfileFolder()
    net       = DataTypes.toNetwork(stationName)
    sta       = DataTypes.toStation(stationName)
    keyfile   = KeyFileObj         (keyfolder, net,sta)
    sta       = keyfile.read()

    if sta == None : site = None
    else :           site = sta.site + '(' + sta.provider + ')'

    return site


def isIRIS(dirName=None, net=None, station=None, fullName=None) :

    if dirName == None : dir = Globals.KeyfileFolder()
    else :               dir = dirName

    provider = getProvider(dir, net,station, fullName)

    if provider == PROV_IRIS   : return True
    if provider == PROV_WEB_DC : return False

    return False # ???
    #Debug.assert1(False, 'Invalid provider in keyfile ' + file.fullName, 'Perhaps old program version') ???
    xxx
# --------------------------------------------------------------------------------------------------

def getIrisMask(dirName, stations) :

    mask = []

    for station in stations :
        #print station
        b = isIRIS(dirName, fullName = station)
        mask.append(b)

    return mask

# --------------------------------------------------------------------------------------------------

def checkVersion(dirName, net=None, station=None, fullName=None) :

    prov = getProvider(dirName, net,station, fullName)

    if prov == None :
       return Logfile.error('Invalid provider in keyfile ' + fullName)

    if prov == PROV_IRIS or prov == PROV_WEB_DC : return True

    return Logfile.error('Invalid provider <' + prov + '> in keyfile station_' + fullName,
                          'Perhaps old program version')
