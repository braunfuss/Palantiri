#!/usr/bin/python -W ignore::DeprecationWarning

import os
import sys


# add local directories to import path

sys.path.append('../tools/')
sys.path.append('../Common/')

from   optparse import OptionParser
import fnmatch
import urllib
import locale
import logging
import multiprocessing
import time
import getpass
import signal
import itertools

from   ConfigParser import SafeConfigParser
from   lxml         import etree

import obspy.core
from obspy.clients.arclink    import Client

import cPickle as pickle

#      Import from Common

import Globals                     # Own global data
import Basic                       # Own module with basic functions
import Logfile                     # Implements logfile
import ConfigFile                  # Semantic of config file entries

import DataTypes
from   DataTypes import Station

#      Import from Waveform

import KeyFile                     # Class KeyFile
import DataDir
import Server                      # Routines for server implementation
import WebDC
from   Version  import  VERSION_STRING

options = None

# -------------------------------------------------------------------------------------------------

g_Diff        = -1
g_MetaInfo    = []

ARCLINK_META   = 'ARCLINK META-DATA'
IRIS_META      = 'IRIS META-DATA'

def metaFileName(options) :

    eventpath = options.path
    metafile  =('metainfo-%s-%s.meta') %(options.year, options.day)

    return os.path.join(eventpath, metafile)

# -------------------------------------------------------------------------------------------------

def parseMetaInfoFile(metafilepath):

    Logfile.add(' ', 'Parsing MetaInfoFile', metafilepath)
    MetaL = []

    try:
        fobj = open(metafilepath,'r')

        for i in fobj:
            line = i.split()
            net  = line[0]; sta  = line[1]; loc  = line[2]; comp = line[3]
            lat  = line[4]; lon  = line[5]; ele  = line[6]; dip  = line[7]
            azi  = line[8]; gain = line[9]; inst = line[10]

            MetaL.append(Station(net,sta,loc,comp,lat,lon,ele,dip,azi,gain,inst))

        Logfile.add(str(len(MetaL)) + ' ENTRIES IN METAFILE FOUND')

    except:
        Logfile.add('METAFILE NOT READABLE', metafilepath)

    return MetaL

# -------------------------------------------------------------------------------------------------

def parseSDS(p):

    Logfile.add(' ')
    Logfile.red('Parsing SDS STRUCTURE ')

    L = []
    p = os.path.join(p,'data')

    for root,dirs,files in os.walk(p):
        for i in files:
            fullName = os.path.join(root,i)                 #hs+

            if os.stat(fullName).st_size == 4 :
               continue                                      #hs- : skip unused channel

            if fnmatch.fnmatch(i, '*.' + options.day):
                    line = str.split(i,'.')
                    net  = line[0]
                    sta  = line[1]
                    loc  = line[2]
                    comp = line[3]
                    L.append(Station(net,sta,loc,comp))

    Logfile.red(str(len(L)) + ' ENTRIES IN SDS STRUCTURE FOUND')
    return L

# -------------------------------------------------------------------------------------------------

def cmpMeta(MetaList, SDSList):
    global g_Diff

    SList = []
    MList = []
    Miss  = []

    for i in SDSList:
        snameSDS = i.net+'.'+i.sta+'.'+i.loc+'.'+i.comp
        SList.append(snameSDS)

    for j in MetaList :
        if j.loc == '--' : j.loc = ''

        snameMeta = j.net+'.'+j.sta+'.'+j.loc+'.'+j.comp
        MList.append(snameMeta)

    t = [item for item in SList if not item in MList]

    g_Diff = len(SList) - len(MList)

    Logfile.add(str(len(SList)+1) + ' ENTRIES IN SDSLIST')               #hs
    Logfile.add(str(len(MList))   + ' ENTRIES IN METALIST')
    Logfile.add(str(len(t))       + ' ENTRIES IN SEARCHDIFF')
    Logfile.add(str(g_Diff)       + ' CALC DIFF')

    for a in t:
        line = a.split('.')
        net  = line[0]
        sta  = line[1]
        loc  = line[2]
        comp = line[3]
        Miss.append(Station(net,sta,loc,comp))

    return Miss

# -------------------------------------------------------------------------------------------------

def makeTime(sdsyear,sdsday):

    t_begin = obspy.core.utcdatetime.UTCDateTime(year=int(sdsyear), julday=int(sdsday),iso8601=True)
    t_end = t_begin + 23 * 60

    iris_begin = t_begin.format_iris_web_service()
    iris_end   = t_end.format_iris_web_service()

    geofon_begin = t_begin.format_arclink()
    geofon_end = t_end.format_arclink()

    time_d = {'g_begin':geofon_begin,'g_end':geofon_end,'i_begin':iris_begin,'i_end':iris_end,'go_begin':t_begin,'go_end':t_end}
    return time_d

# -------------------------------------------------------------------------------------------------

def getMetaInfo(filename):

    Logfile.red('PARSING METAINFOFILE %s' %(filename))
    locale.setlocale(locale.LC_ALL, 'C')
    meta_dict = {}

    sp  = Parser(filename)
    d   = sp.stations
    t   = d[0]
    gain= locale.format('%.2f',t[3].sensitivity_gain,1)
    lat = locale.format('%.4f',t[1].latitude,1)
    lon = locale.format('%.4f',t[1].longitude,1)
    ele = locale.format('%.2f',t[1].elevation,1)
    dip = locale.format('%.2f',t[1].dip,1)
    azi = locale.format('%.2f',t[1].azimuth,1)

    if gain == 0 : meta_dict = {}
    else:          meta_dict = {'lat':lat,'lon':lon,'ele':ele,'dip':dip,'azi':azi,'gain':gain}

    return meta_dict

# -------------------------------------------------------------------------------------------------
#hs : ersetzt durch getArclinkInst_2()

def getArclinkInst(station,usermail,pwdkeys):

    inst = 'no_instrument_available'

    #if station.loc == '' : station.loc='--'

    sn = station.net+'.'+station.sta+'.'+station.loc+'.'+station.comp
    print 'SN ',sn

    client = Client(user=usermail,dcid_keys=pwdkeys)
    inv    = client.getInventory(station.net, station.sta, station.loc, station.comp, instruments=True)
    print 'INV: ',inv

    try:
        stats = inv[sn]
        print 'STATS ',stats

        t    = stats[0]
        z    = t['paz']
        inst = z['name']
        inst = inst.replace(' ','_')

    except:
        inst = 'no_instrument_available'

    return inst

# -------------------------------------------------------------------------------------------------
#hs+

def getFromKeyFile(sta) :

    file = KeyFile.KeyFileObj(dirName=None, net=sta.net, station=sta.sta)
    sta1 = file.read()

    if sta1 == None :
       print 'No keyfile'
       return station

    sta2     = sta
    sta2.lon = sta1.lon
    sta2.lat = sta1.lat
    sta2.ele = sta1.ele

    return sta2


def getArclinkInst_2(station, usermail, pwdkeys):

    #if station.loc == '--' : loc = ''

    sta2      = station
    sta2.lon  = -1
    sta2.lat  = -1
    sta2.ele  = -1
    sta2.gain = -1
    sta2.dip  = -1
    sta2.azi  = -1
    sta2.inst = 'no_instrument_available'

    loc = station.loc

    if loc == '00' or loc == '01' or loc == '10' or loc == '11' :     # no entry for this loc
       print 'Use key file'
       return getFromKeyFile(sta2)

    loc    = ''

    client = Client(user=usermail, dcid_keys=pwdkeys)

    try    : inv = client.getInventory(station.net, station.sta, loc, station.comp, instruments=True)
    except : return getFromKeyFile(sta2)

    #inv = WebDC.getInventory(station, usermail, pwdkeys)  # ??? noch einbauen

    print 'INV: ',inv
    keys = inv.keys()

    if len(keys) == 0 :
       print 'No keys'
       return None

    if Globals.isDebug :               #    save inventory

       file    = Server.InventoryFileName(station, '.txt')
       fp      = open(file, 'w')
       lines   = []

       for k in keys :
          lines.append('key = ' + str(k) + '\n')
          lines.append('inv = ' + str(inv[k]) + '\n')

       Basic.writeTextFile(file, lines)
    #endif

    inv1     = inv  [station.net + '.' + station.sta]
    sta2.lon = inv1 ['longitude']
    sta2.lat = inv1 ['latitude']
    sta2.ele = inv1 ['elevation']

    key  = station.net + '.' + station.sta + '.' + loc + '.' + station.comp
    inv1 = inv [key]

    z  = inv1 [0]
    z1 = z ['paz']

    if 'gain' in z1 : sta2.gain = z1 ['gain']

    if 'name' in z1 :
       inst      = z1 ['name']
       sta2.inst = inst.replace(' ','_')

    sta3 = sta2
    return sta3

#hs-
# -------------------------------------------------------------------------------------------------

def arclinkRequest_sc(station,begin,end):

    sname = station.net+'.'+station.sta+'.'+station.loc+'.'+station.comp
    print 'STARTING ARCLINK REQUEST FOR STATION ' + sname

    dt = obspy.core.utcdatetime.UTCDateTime(begin)
    et = obspy.core.utcdatetime.UTCDateTime(end)

    reqfile = 'ARC-'+sname+'-'+station.loc+'.txt'
    fname   = os.path.join(os.getcwd(),reqfile)
    fobj    = open(fname,'w')

    if station.loc == '--' : station.loc = ''

    content = str(dt.format_arclink())+' '+str(et.format_arclink())+' '+station.net+' '+station.sta+' '+station.comp+' '+station.loc+'\n'
    Logfile.add('%s %s %s' %(sname, str(dt.format_arclink()) ,str(et.format_arclink())))
    fobj.write(content)
    fobj.close()

    output = sname+'-ARCLINK.dataless'
    cmd    = ''
    email  = 'ehlert@geo.uni-potsdam.de'

    if fnmatch.fnmatch(station.net, 'MN'):
        cmd = 'arclink_fetch -qqq -p -a webdc.eu:18002 -w dcidpasswords.txt -k dseed -u ' + email + ' -o '+output+' '+fname
    else:
        cmd = 'arclink_fetch -qqq -a webdc.eu:18002    -w dcidpasswords.txt -k dseed -u ' + email + ' -o '+output+' '+fname

    Logfile.add('%s '%(cmd))
    os.system(cmd)
    info = ''

    try:
        mdict = getMetaInfo(output)
        inst  = getArclinkInst(station)

        if station.loc == '' : station.loc = '--'

        info = '{0:3}{1:6}{2:3}{3:4}{4:12}{5:12}{6:10}{7:10}{8:10}{9:20}{10:50}'.format(station.net.strip(),station.sta.strip(), \
        station.loc.strip(), station.comp.strip(), \
        mdict['lat'], mdict['lon'], mdict['ele'], mdict['dip'],  mdict['azi'], mdict['gain'], inst)

        #print 'ARCLINK METAINFO FOR STATION %s' %(info)
        os.remove(output)
        os.remove(reqfile)

    except:
        print 'NO ARCLINK METAINFO FOR STATION'

    try:
        os.remove(reqfile)
        os.remove(output)
    except:
        ''

    return info

# -------------------------------------------------------------------------------------------------

def arclink_Request_obspy(station, begin, end, usermail, pwdkeys):

    sname = station.net+'.'+station.sta+'.'+station.loc+'.'+station.comp
    print 'STARTING ARCLINK REQUEST FOR STATION ' + sname

    output = sname + '-ARCLINK.dataless'
    client = Client(user=usermail,dcid_keys=pwdkeys,debug=False)
    info   = ''

    if station.net == 'MN':
        info = arclinkRequest_sc(station,begin,end)    #hs : Kommt wohl nicht vor - wird immer von Iris geholt

    else:
        try:
#           client.saveResponse(output, station.net, station.sta, station.loc, station.comp, begin, end) #hs

#           mdict = getMetaInfo(output)                                  #hs
#           inst  = getArclinkInst(station, usermail, pwdkeys)         #hs
            sta2  = getArclinkInst_2(station, usermail, pwdkeys)         #hs
            #sta2.print1()

            if station.loc == '' : station.loc = '--'

            form = '{0:3}{1:6}{2:3}{3:4}{4:12}{5:12}{6:10}{7:10}{8:10}{9:20}{10:50}'

            info = form.format(station.net.strip(), station.sta.strip(), station.loc.strip(),
                                station.comp.strip(),
#           mdict['lat'], mdict['lon'], mdict['ele'], mdict['dip'], mdict['azi'], mdict['gain'],inst)      #hs
            str(sta2.lat),  str(sta2.lon), str(sta2.ele), str(sta2.dip), str(sta2.azi),               #hs \
            str(sta2.gain), str(sta2.inst))                                                              #hs

            #logfile.red('ARCLINK METAINFO FOR STATION %s' %(info))
#           os.remove(output)                        #hs
            print ARCLINK_META                          #hs
            return info

        except:
            print 'Exception 1 in arclink_Request_obspy'
            info = ''
            return info
    try:
        return info
        os.remove(output)

    except:
        print 'Exception 2 in arclink_Request_obspy'
        pass

    print 'len=', len(info)
    return info

# -------------------------------------------------------------------------------------------------
#hs+

def xmlValue(token, line) :
    s1 = line.split('>')
    s2 = s1[1].split('<')
    return s2[0]


def xmlString(token, line, default='???') :

    if not token in line : return default
    else :                 return xmlValue(token, line)

def xmlFloat(token, line, default=-1) :

    if not token in line : return default
    else :                 return float(xmlValue(token, line))

# -------------------------------------------------------------------------------------------------

def extractIrisData(station, xmlFileData) :

    if Globals.isDebug :               #    save xml file

       file    = Server.xmlFileName(station, '.xml')
       fp      = open(file, 'w')
       fp.write(xmlFileData)
       fp.close()
    #endif

    #    extract data
    #
    lines = xmlFileData.split('\n')

    lat  = -1; lon = -1; ele = -1; dip = -1; azi = -1; gain = -1
    inst = '???'

    for i in range(len(lines)) :
       s = lines[i]

       if '<Latitude>'              in s : lat  = xmlFloat('<Latitude>',  s)
       if '<Longitude>'             in s : lon  = xmlFloat('<Longitude>', s)
       if '<Elevation>'             in s : ele  = xmlFloat('<Elevation>', s)
       if '<Azimuth>'               in s : azi  = xmlFloat('<Azimuth>',   s)
       if '<Dip>'                   in s : dip  = xmlFloat('<Dip>',       s)
       if '<InstrumentSensitivity>' in s : gain = xmlFloat('<Value>', lines[i+1])
       if '<Sensor>'                in s : inst = xmlString('<Type>',  lines[i+1])
    #endfor

    sta = Station(station.net, station.sta, station.loc, station.comp,lat,lon,ele,dip,azi,gain,inst)
    return sta

# -------------------------------------------------------------------------------------------------
# Infos siehe : http://service.iris.edu/fdsnws/station/1/

def irisParams(station, begin, end) :

    params =  'net='  + station.net
    params += '&sta=' + station.sta
    params += '&loc=' + station.loc
    params += '&cha=' + station.comp
    params += '&starttime=' + begin + '&endtime=' + end
    params += '&level=resp'
    return params

#hs-

# -------------------------------------------------------------------------------------------------

def irisRequest(station,begin,end):

    if station.loc == '' : station.loc='--'
#hs+
#   url = 'http://www.iris.edu/ws/station/query?net=' + \
#   station.net+'&sta='+station.sta+'&loc='+station.loc+'&cha='+station.comp+'&starttime='+begin+'&endtime='+end + '&level=resp'

    url = 'http://service.iris.edu/fdsnws/station/1/query?' + irisParams(station,begin,end)
#hs-

    meta = ''

    try:
#       f    = urllib.urlopen(url)               #hs+
#       doc  = etree.parse(f)
#       root = doc.getroot()
#
#       lat=''
#       lon=''
#       ele=''
#       dip=''
#       azi=''
#       gain=''
#       inst=''
#
#       for child in root:
#          for a in child:
#              for b in a:
#                for c in b:
#                   for d in c:
#                       for e in d:
#                          if e.tag == '{http://www.data.scec.org/xml/station/}Lat':         lat = e.text
#                          if e.tag == '{http://www.data.scec.org/xml/station/}Lon':         lon = e.text
#                          if e.tag == '{http://www.data.scec.org/xml/station/}Elevation':   ele = e.text
#                          if e.tag == '{http://www.data.scec.org/xml/station/}Azimuth':     azi = e.text
#                          if e.tag == '{http://www.data.scec.org/xml/station/}Dip':         dip = e.text
#
#                          for f in e:
#                             if f.tag == '{http://www.data.scec.org/xml/station/}SensitivityValue' : gain = f.text
#                             if f.tag == '{http://www.data.scec.org/xml/station/}EquipType' :        inst = f.text

        xmlFileData = urllib.urlopen(url).read()

        if len(xmlFileData) == 0 :
           return ''

        sta2 = extractIrisData(station, xmlFileData)         # sta2 is of type Station

        lat  = str(sta2.lat); lon = str(sta2.lon);  ele  = str(sta2.ele)
        dip  = str(sta2.dip); azi = str(sta2.azi);  gain = str(sta2.gain)
        inst = str(sta2.inst)                                                   #hs-

        inst = inst.replace(' ','_')

        if gain != '' :
           gain = float(gain)
           locale.setlocale(locale.LC_ALL,'C')
           x = locale.format('%.2f',gain,1)

        meta = '{0:3}{1:6}{2:3}{3:4}{4:12}{5:12}{6:10}{7:10}{8:10}{9:20}{10:50}'.format(station.net.strip(), \
        station.sta.strip(),station.loc.strip(),station.comp.strip(),lat.strip(),lon.strip(),ele.strip(),     \
        dip.strip(),azi.strip(),x.strip(),inst.strip())

        print IRIS_META

    except:
        print 'Exception'

    return meta

# -------------------------------------------------------------------------------------------------

def makeRequest_old(i, time_d, flag, usermail, pwdkeys):

        sname = i.net+'.'+i.sta+'.'+i.loc+'.'+i.comp
        meta  = irisRequest(i, time_d['i_begin'], time_d['i_end'])

        if meta == '' :
           meta = arclink_Request_obspy(i,time_d['go_begin'],time_d['go_end'], usermail, pwdkeys)

        serializeMeta(meta,flag)

        if meta != '' :
           print '#meta = '                              #hs
           print meta                                    #hs


def makeRequest(i, time_d, flag, usermail, pwdkeys):

    sname  = i.net+'.'+i.sta+'.'+i.loc+'.'+i.comp
    p      = KeyFile.getProvider(net=i.net, station=i.sta)
    meta   = ''
    t1     = time_d ['i_begin']
    t2     = time_d ['i_end']

    try :
       if   p == KeyFile.PROV_IRIS :    meta = irisRequest(i, t1, t2)

       elif p == KeyFile.PROV_WEB_DC :
          try :    meta = arclink_Request_obspy(i, t1, t2, usermail, pwdkeys)
          except : meta = ''

          if meta == '' :
             try :    meta = irisRequest(i, t1, t2)
             except : meta = ''

       else :   # Keyfile not found
          try :    meta = irisRequest(i, t1, t2)
          except : meta = ''

          if meta == '' :
             try :    meta = arclink_Request_obspy(i, t1, t2, usermail, pwdkeys)
             except : meta = ''

    except : meta = ''

    if meta != '' :
       print '#meta = '
       print meta

# -------------------------------------------------------------------------------------------------

def writeList(L,PreList,numproc,metafilepath):

    fname = metafilepath
    Logfile.red('WRITING METADATAFILE ----------->  %s' %(fname))
    print options.path

    fobj  = open(os.path.join(os.getcwd(),fname),'w')
    WLIST = []

    for j in PreList:
        if j.loc == '' : j.loc = '--'

        form = '{0:3}{1:6}{2:3}{3:4}{4:12}{5:12}{6:10}{7:10}{8:10}{9:20}{10:50}'
        meta = form.format(j.net,j.sta,j.loc,j.comp,j.lat,j.lon,j.ele,j.dip,j.azi,j.gain,j.inst)
        WLIST.append(meta)

    for i in L : WLIST.append(i)

    WLIST = sorted(WLIST)

    fobj.write('\n'.join(WLIST))
    fobj.close()

# -------------------------------------------------------------------------------------------------

def globalConf():

    cDict  = {}
    parser = SafeConfigParser()
    parser.read('../global.conf')

    for section_name in parser.sections() :
        for name, value in parser.items(section_name) :
            cDict[name]=value

    return cDict

# -------------------------------------------------------------------------------------------------

g_LastStation = None

def checkProcessError(stationName, nErrors, lines, execTime) :   # ??? execTime

    global g_MetaInfo, g_LastStation

    errCode = Server.HAS_NO_DATA
    errCode = Server.RETRY_IT

    for lineNr in range(len(lines)) :
       line  = lines [lineNr]
       isEnd = False
       s     = None

       # UserWarning: MAX_REQUESTS exceeded - breaking current request loop

       if 'MAX_REQUESTS' in line :
          errCode = Server.RETRY_IT

          s =  'UserWarning: MAX_REQUESTS exceeded - breaking current request loop'
          s += '(' + str(nErrors) + ')'
          #isEnd = True

       elif 'deprecated' in line : s = ' '             # ignore ObsPyDeprecation Warning   #15.7.2016

       elif Logfile.MSG_TOKEN in line :  s = line
       elif 'python'  in line :          s = line

       elif ARCLINK_META in line or IRIS_META in line :
          name = DataTypes.toNetAndStation(stationName)

          if g_LastStation == None or g_LastStation != name :
             g_LastStation = name
             s = KeyFile.getSite(stationName)

       elif Server.CLIENT_ABORT_MSG in line :
          errCode = Server.RETRY_IT
          s       = line

       elif '#meta' in line  :                      # station has data
          errCode = Server.HAS_DATA
          isEnd   = True

          s = lines [lineNr+1]
          g_MetaInfo.append(s)

       elif 'Traceback' in line :
          sn = []

          for i in range(0,300) :
              if lineNr+i >= len(lines) : break         #10.12.2015
              if 'KeyboardInterrupt' in lines [lineNr+i] : sn = []; break
             #if lineNr+i >= len(lines) : break         #10.12.2015

              sn.append(lines [lineNr+i])
          #endfor

          if Server.checkIsTimeOut(stationName, sn) :       # Traceback shows timeout
             Logfile.error('Retry access later')
             errCode = Server.RETRY_IT

          else :                                         # Traceback --> log
             Server.printLines(stationName, sn, onlyErrorLog=True)

          isEnd = True
       #endif

       if s != None : Server.printMsg(stationName, s)
       if isEnd     : break
    #endwhile

    return errCode

# -------------------------------------------------------------------------------------------------

def writeMetaInfo(metaInfo) :

    file     = metaFileName(options)
    oldLines = Basic.readTextFile(file)
    info2    = sorted(itertools.chain(oldLines, metaInfo))
    lines    = []

    for line in info2  : lines.append(line)

    Basic.writeTextFile(file, lines)
    Logfile.add('Write Meta information to file ', file, ' ')

# -------------------------------------------------------------------------------------------------

def buildStationList(options) :

    metafile  =('metainfo-%s-%s.meta') %(options.year, options.day)
    metaf     = os.path.join(options.path, metafile)

    time_d = makeTime(options.year, options.day)
    SDS    = parseSDS(options.path)

    #   Select used channels
    #
    names = []

    for s in SDS : names.append(s.fullName())

    old = False

    if not old : SDS_2 = SDS
    else :
       selection = ['BHE','BHN','BHZ', 'BH1','BH2','BHZ', 'HHE','HHN','HHZ']    # use only this
       mask      = Basic.stringsEndsWith(names, selection)

       SDS_2 = []

       for i in range(len(mask)) :
           if mask[i] : SDS_2.append(SDS[i])
    #endif

    #
    #
    META   = parseMetaInfoFile(metaf)
    D      = cmpMeta(META, SDS_2)

    waiting = []
    n       = len(D)

    #if Globals.isDebug :
    #   if n > 200 : n = 200

    for i in range(n) :  waiting.append(D[i].fullName())

    return waiting

# -------------------------------------------------------------------------------------------------

def init(isClient) :

    Globals.isClient = isClient

    if not isClient :
       if not Logfile.init(startMsg = VERSION_STRING) : return False

    return Globals.init()

# --------------------------------------------------------------------------------------------------

def checkConfigFile(conf) :

    mail          = ConfigFile.mail
    pwd           = ConfigFile.pwd
    keyfilefolder = ConfigFile.keyfilefolder

    keyList = [mail, pwd, keyfilefolder]
    return ConfigFile.checkKeys(conf, keyList)

# -------------------------------------------------------------------------------------------------

SERVER_NAME = 'meta'

def startIrisServer(stations, args) :

    ctrl = Server.ServerCtrl(nRetries = 4, nParallel=10, waitTime=0.1, printStat=False)
    srv  = Server.ServerBase(SERVER_NAME, checkProcessError, ctrl)
   #if WINDOWS : srv.control.ClientProc = MainProc

    return srv.run(stations, args)


def startGeofonServer(stations, args) :

    ctrl = Server.ServerCtrl(nRetries = 4, nParallel=4, waitTime=0.5, printStat=False)
    srv  = Server.ServerBase(SERVER_NAME, checkProcessError, ctrl)
   #if WINDOWS : srv.control.ClientProc = MainProc

    return srv.run(stations, args)

# -------------------------------------------------------------------------------------------------

def startServer(stationList, options) :

     network    = options.network
     mask       = KeyFile.getIrisMask(None, stations=stationList)
     irisList   = Basic.selectStrings(stationList, mask)
     geofonList = Basic.selectStrings(stationList, Basic.Not(mask))

     Conf = Globals.ConfigDict
     args = Server.joinClientArgs([Conf ['pwd'], Conf ['mail']])

     if not network or network == 'iris' :
        if len(irisList) == 0 :
           Logfile.add('All iris entries set')

        else :
           if not startIrisServer(irisList, args) :
              return True                              # aborted with ctrl c

     if not network or network == 'geofon' :
        if len(irisList) == 0 :
           Logfile.add('All geofon entries set')

        else :
           startGeofonServer(geofonList, args)
           return True

     if network and network != 'iris' and network != 'geofon' :
        if not DataDir.isNetwork(network) :
           return Logfile.error('Illegal network name <' + network + '>')

        list2 = DataTypes.selectNetwork(irisList, network)

        if len(list2) > 0 :
           startIrisServer(list2, args)
           return True

        list2 = DataTypes.selectNetwork(geofonList, network)

        if len(list2) > 0 :
           startGeofonServer(list2, args)
           return True

        Logfile.add('All network enties set')
     #endif

     return False       # nothing done

# --------------------------------------------------------------------------------------------------
#
#   Client routine
#

class MetaDataClient(Server.ClientBase) :

    def __init__(self, options) :

        self.options = options
        self.args    = Server.splitClientArgs(options.station)
        self.station = self.args[0]

        Server.ClientBase.__init__(self, self.station)

    def _run(self) :          # called from ClientBase

       #Conf    = globalConf()
       #pwd     = Conf ['pwd']
       #mail    = Conf ['mail']
        pwd     = self.args[1]
        mail    = self.args[2]
        pwdDict = {}

        for i in pwd.split(','):
            pwdDict [i.split(':')[0]] = i.split(':')[1]

        time_d  = makeTime(self.options.year, self.options.day)
        names   = self.station.split('.')
        station = Station(names[0], names[1], names[2], names[3])

        makeRequest(station, time_d, 0, mail, pwdDict)

#endclass MetaDataClient

# --------------------------------------------------------------------------------------------------

def run_parallel(options) :

    if options.station :                                     # Client part
       Globals.setEventDir(options.path)

       if not init(True) : return False

       clt = MetaDataClient(options)
       clt.run()

    else :                                                    # Server part
       Basic.checkExistsDir(options.path, isAbort=True)
       Globals.setEventDir(options.path)

       if not init(False) : return False

       checkConfigFile(Globals.ConfigDict)

       #   Build station list
       stationList = sorted(buildStationList(options))

       if len(stationList) == 0 :
          if g_Diff != 0 : msg = 'No stations with data found : Run getdata before'
          else :           msg = 'MetaInfoFile complete'

       else :
          #   Run server
          #
          if not startServer(stationList, options) :
             msg = 'Program finished'

          else :
             writeMetaInfo(g_MetaInfo)
             buildStationList(options)

             if g_Diff != 0 : msg = 'Program finished'
             else :           msg = 'Program finished - All stations found'
       #endif

       Logfile.showLabel(msg)
    #endif server

# -------------------------------------------------------------------------------------------------
#  Aufruf : python  tools/ev_meta_mt4.py  -f  <event dir>
#
def MainProc() :

    global options

    p = OptionParser(usage="%prog -p PATH -y YEAR -d JULIANDAY")

    p.add_option("-p","--path",   type="string", dest="path",    help="time")
    p.add_option("-y","--year",   type="string", dest="year",    help="year")
    p.add_option("-d","--day",    type="string", dest="day",     help="day")
    p.add_option("-x","--dummy",  type="string", dest="station", help="dummy")   #hs : client flag
    p.add_option("-n","--dummy2", type="string", dest="network", help="dummy2")  #hs : use single network

 (options, args) = p.parse_args()
    run_parallel(options)

if __name__ == "__main__":

    MainProc()
