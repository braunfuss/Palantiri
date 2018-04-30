import re
import os
import fileinput
from obspy.taup.taup import getTravelTimes
from obspy.taup.taup import locations2degrees
from obspy.core.utcdatetime import UTCDateTime
from obspy.arclink.client import Client
from obspy.core.stream import Stream
from obspy.core import read
from obspy.signal.filter import bandpass as bp
from obspy.signal.filter import highpass as hp
from obspy.signal.filter import lowpass as lp
import logging
import fnmatch
import numpy as np
import time
from optparse import OptionParser
import shutil
import math

logger = logging.getLogger('AT')
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

parser = OptionParser(usage="%prog -f Eventpath ")
parser.add_option("-f", "--evpath", type="string", dest="evpath", help="evpath")

(options, args) = parser.parse_args()

if options.evpath == None:
        parser.error("non existing eventpath")

evpath = options.evpath



class Station(object):
    def __init__(self, net, sta, loc, comp,lat=0,lon=0,ele=0,dip=0,azi=0,gain=0):
        self.net = net
        self.sta = sta
        self.loc = loc
        self.comp = comp
        self.lat = lat
        self.lon = lon
        self.ele = ele
        self.dip = dip
        self.azi = azi
        self.gain = gain
    def getName(self):
        return self.net+'.'+self.sta+'.'+self.loc+'.'+self.comp

class GridElem (object):
    def __init__(self, lat, lon, depth,tt,delta):
        self.lat = lat
        self.lon = lon
        self.depth = depth
        self.tt = tt
        self.delta = delta
        
class TTTGrid(object):
    def __init__(self, dimZ, mint,maxt,Latul,Lonul,Lator,Lonor,GridArray):
        self.dimZ = dimZ
        self.mint = mint
        self.maxt = maxt
        self.Latul = Latul
        self.Lonul = Lonul
        self.Lator = Lator
        self.Lonor = Lonor
        self.GridArray = GridArray

def readEvent(EventPath):

    Origin = {}
    logger.info('\033[31m Parsing EventFile \033[0m \n')
    try:
        for i in os.listdir(EventPath):
            if fnmatch.fnmatch(i, '*.origin'):
                evfile = os.path.join(EventPath,i)
        fobj = open(evfile, 'r')
        
        for line in fobj:
            line = str.split(line, '=')
            key  = (line[0].replace('\'', '')).strip()
            value = line[1].replace('\n','').strip()

            Origin[key] = value
    except:
        logger.info('\033[31m File Error Exception \033[0m \n') 

    return Origin

def readConfig(EventPath):
    
    Config = {}
    logger.info('\033[31m Parsing ConfigFile \033[0m \n')
    try:
        for i in os.listdir(EventPath):
            if fnmatch.fnmatch(i, '*.config'):
                evfile = os.path.join(EventPath,i)
        fobj = open(evfile, 'r')
        
        for line in fobj:
            if line[0] != '#' and len(line) > 1:
                line = str.split(line, '=')
                key  = (line[0].replace('\'', '')).strip()
                value = line[1].replace('\n','').strip()

                Config[key] = value

    except:
        logger.info('\033[31m File Error Exception \033[0m \n')

    return Config

def readMetaInfoFile(EventPath):

    logger.info('\033[31m Parsing MetaInfoFile \033[0m \n')
    try:
        for i in os.listdir(EventPath):
            if fnmatch.fnmatch(i, '*.meta'):
                evfile = os.path.join(EventPath,i)

        MetaL = []
    
        fobj = open(evfile,'r')
        for i in fobj:
            line = i.split()
            net = line[0]
            sta = line[1]
            loc = line[2]
            comp = line[3]
            lat = line[4]
            lon = line[5]
            ele = line[6]
            dip = line[7]
            azi = line[8]
            gain = line[9]
            
            if fnmatch.fnmatch(comp,'*HZ'):
                MetaL.append(Station(net,sta,loc,comp,lat,lon,ele,dip,azi,gain))

        logger.info('\033[31m %d ENTRIES IN METAFILE FOUND \033[0m \n' % (len(MetaL)))
    except:
        logger.info('\033[31m METAFILE NOT READABLE \033[0m \n')

    return MetaL

def createFolder(EventPath):
    
    Folder = {}
    logger.info('\033[31m Create working environment \033[0m \n')
    if os.access(os.getcwd(), os.W_OK):
        
        basedir = os.path.join(EventPath,'work')
        sembdir = os.path.join(basedir, 'semblance')
        ascdir  = os.path.join(basedir, 'asc')
        mseeddir  = os.path.join(basedir, 'mseed')
        
        Folder['base'] = basedir
        Folder['semb'] = sembdir
        Folder['asc']  = ascdir
        Folder['mseed'] = mseeddir
        
        for key in Folder:
            if os.access(Folder[key],os.F_OK) == False:
                os.makedirs(Folder[key])
                
    else:
        print"no write permissions"
        
    Folder['config'] = os.path.join(os.getcwd(),'skeleton')
    return Folder

def cpSkeleton(FolderDict,ConfigDict):

    logger.info('\033[31m Copy needed processing scripts \033[0m \n')
    if os.access(os.getcwd(), os.W_OK):
        try:
            shutil.copyfile(os.path.join(FolderDict['config'],ConfigDict['plot1']), os.path.join(FolderDict['semb'],ConfigDict['plot1']))
            shutil.copyfile(os.path.join(FolderDict['config'],ConfigDict['plot2']), os.path.join(FolderDict['semb'],ConfigDict['plot2']))
            shutil.copyfile(os.path.join(FolderDict['config'],ConfigDict['rerun']), os.path.join(FolderDict['semb'],ConfigDict['rerun']))
        except:
            logger.info('\033[31m Copy processing scripts Error \033[0m \n')



def filterStations(StationList,Config,Origin,network):
    
    F = []
    minDist = int(Config['minDist'])
    maxDist = int(Config['maxDist'])
    
    o_lat = float(Origin['lat'])
    o_lon = float(Origin['lon'])
       
    logger.info('\033[31m Filter stations with configured parameters \033[0m')
    
    for i in StationList:
        if i.loc == '--':
            i.loc=''
        streamID = i.net+'.'+i.sta+'.'+i.loc+'.'+i.comp
       

        for j in network:

            if fnmatch.fnmatch(streamID, j):
                sdelta = locations2degrees(o_lat, o_lon, float(i.lat), float(i.lon))

                if sdelta > minDist and sdelta < maxDist:
                    F.append(Station(i.net,i.sta,i.loc,i.comp,i.lat,i.lon,i.ele,i.dip,i.azi,i.gain))


    logger.info('\033[31m %d STATIONS LEFT IN LIST \033[0m'% len(F))
    return F


def calcTTT(Config,StationList,Origin):
    
    dimX = int(Config['dimX']) 
    dimY = int(Config['dimY']) 
    gridspacing = float(Config['gridspacing'])
    
    o_lat = float(Origin['lat'])
    o_lon = float(Origin['lon'])
    o_depth = float(Origin['depth'])

    logger.info(' BOUNDING BOX DIMX: %d  DIMY: %d  GRIDSPACING: %f \n'%(dimX,dimY,gridspacing))    
    
    
    oLator = o_lat + dimX/2
    oLonor = o_lon + dimY/2
    oLatul = 0
    oLonul = 0
    mint= 100000
    maxt=-100000
    
    TTTGridMap = {}
    
    for station in StationList:
        
        GridArray = {}
        streamID = station.net+'.'+station.sta+'.'+station.loc+'.'+station.comp
        sdelta = locations2degrees(float(o_lat), float(o_lon), float(station.lat), float(station.lon))
              
        logger.info(' STATION: %s --> DELTA: %f'% (streamID,sdelta))
        
        z=0
        for i in xrange(dimX):
            
            oLatul = o_lat -((dimX-1)/2)*gridspacing + i*gridspacing
            if z == 0 and i == 0:
                Latul = oLatul
            o=0    
            for j in xrange (dimY):
                
                oLonul = o_lon -((dimY-1)/2)*gridspacing + j*gridspacing
                if o==0 and j==0:
                    Lonul = oLonul
                
                de = locations2degrees(float(oLatul), float(oLonul), float(station.lat), float(station.lon))
                
                tt = getTravelTimes(delta=de,depth=o_depth,model='ak135')
                #print tt
                if tt[0]['phase_name'] == Config['ttphase']:
                        time = tt[0]['time']
                        #print streamID,oLatul,oLonul,' --------> ' ,time , ' -> \n'
                        
                        GridArray[(i,j)] = GridElem(oLatul, oLonul, o_depth,time,de)
                    
                        if (mint > time):
                                mint = time
                        if (maxt < time):
                                maxt = time        


        TTTGridMap[streamID] = TTTGrid(o_depth,mint,maxt,Latul,Lonul,oLator,oLonor,GridArray)
                                                        
    logger.info('\033[31m MINT: %g  MAXT: %f \033[0m'% (mint,maxt))

    return mint, maxt,TTTGridMap


def calculateTimeWindows(mint,maxt,Config,Origin):

    winlen = int(Config['winlen'])
    forerun = int(Config['forerun'])
    duration = int(Config['duration'])
    security = int(Config['security'])

    tw = {}

    o_time = UTCDateTime(Origin['time'])
    tw['start'] = UTCDateTime(o_time+(mint-forerun-security))
    tw['end'] = UTCDateTime(o_time+(maxt+duration+winlen+security))
 
    logger.info(' ORIGIN TIME %s'% o_time)
    logger.info(' TIME WINDOW: %s - %s' % (tw['start'], tw['end']) )

    return tw

def readWaveforms(stationList,tw,EventPath,Origin):
    
    time = Origin['time']
    ts = time.split('-')
    year = ts[0].strip()
    month = ts[1]
    ds = ts[2].split(' ')
    day = ds[0]
    
    julday = UTCDateTime(int(year),int(month),int(day)).julday
    
    sdspath = os.path.join(EventPath,'data',year)

    Wdict = {}
    stream = Stream()

    for i in stationList:
        streamID = i.net+'.'+i.sta+'.'+i.loc+'.'+i.comp
        streamData = streamID+'.D.'+str(year)+'.'+str(julday)
        entry = os.path.join(sdspath,i.net,i.sta,i.comp+'.D',streamData)

        st = read(entry,format="MSEED",starttime=tw['start'],endtime=tw['end'],nearest_sample=True)
        
        if len(st.getGaps()) > 0:
            st.merge(method=0, fill_value='interpolate', interpolation_samples=0)
        
        Wdict[streamID] = st
        print st
        stream += st
        logger.info('%s added to StreamList ', (streamID))

    logger.info('%s Streams with available Data ', (len(stationList)))
    
    return stream,Wdict

def writeWaveform(Folder,station,Stream,flag,network):
    
    if flag == 'U':
        name = os.path.join(Folder['mseed'],network+'-'+station+'-U.mseed')
        Stream.write(name,format='MSEED')
        logger.info('unfiltered stream for station %s written '% (station))
        
    elif flag == 'F':
        name = os.path.join(Folder['mseed'],network+'-'+station+'-F.mseed')
        Stream.write(name,format='MSEED')
        logger.info('filtered stream for station %s written '% (station))
        
    elif flag == 'R':
        name = os.path.join(Folder['mseed'],network+'-'+station+'-R.mseed')
        Stream.write(name,format='MSEED')
        logger.info('resampled stream for station %s written '% (station))

def getGain(streamID,MetaDict):
    
    gain=1
    
    for i in MetaDict:
        stream = i.net+'.'+i.sta+'.'+i.loc+'.'+i.comp 
        if fnmatch.fnmatch(streamID, stream):
            gain = i.gain

    return gain



def process(WaveformDict,Config,Folder,network,MetaDict):
    
    new_frequence = float(Config['new_frequence'])
    
    #loop over streammap
    for i in WaveformDict:

            if Config['export_unfiltered'].capitalize() == 'True':
                    writeWaveform(Folder,i,WaveformDict[i],'U',network)

            # TODO implement gain multiplication
            gain = getGain(i,MetaDict)
            print 
            '''
            for j in WaveformDict[i]:
                newtracedata = np.ndarray(shape=(j.stats.npts,1), dtype=float)
                for h in range(j.stats.npts):
                    newtracedata[h] = j.data[h]* 1/float(gain)
                j.data = newtracedata
            '''
            #filtering
            if int(Config['filterswitch']) == 1:
                
                WaveformDict[i][0].filter('bandpass',freqmin = float(Config['flo']),
                           freqmax = float(Config['fhi']),
                           corners = float(Config['ns']),
                           zerophase=bool(Config['zph']))

            elif int(Config['filterswitch']) == 2:
               WaveformDict[i][0].filter("lowpass",
                                      freq=float(Config['l_fc']),
                                      corners=int(Config['l_ns']),
                                      zerophase=bool(Config['l_zph']))

            elif int(Config['filterswitch']) == 3:
                WaveformDict[i][0].filter("highpass", 
                                       freq=float(Config['h_fc']),
                                       corners=int(Config['h_ns']),
                                       zerophase=bool(Config['h_zph']))

            else:
                logger.info('no filter set for station %s '% (i))

            if Config['export_filtered'].capitalize() == 'True':
                    writeWaveform(Folder,i,WaveformDict[i],'F',network)

            #resampling
            logger.info('downsampling %s to %d Hz'% (i,new_frequence))
            factor = WaveformDict[i][0].stats.sampling_rate/new_frequence
            WaveformDict[i][0].decimate(int(factor),strict_length=False, no_filter=False)
            print WaveformDict[i]

            if Config['export_resampled'].capitalize() == 'True':
                writeWaveform(Folder,i,WaveformDict[i],'R',network)

    return WaveformDict

def writeConfig(Config,Origin,Folder):
    
    fobj = open(os.path.join(Folder['semb'],'stations_0.cfg'),'w')
    fobj.write('%source: '+Origin['lat']+' '+Origin['lon']+' '+Origin['depth']+' '+Origin['time']+'\n')
    for i in Config:
        fobj.write('%'+i+': '+Config[i]+'\n')
    fobj.close()

def writeStationFile(StationMetaData,Folder,flag):
    
    name = 'stations_'+str(flag)+'.dat'
    fobj = open(os.path.join(Folder['semb'],name),'w')
    for i in StationMetaData:
        fobj.write('%s %s %s %s\n' %(flag,i.getName(),i.lat,i.lon))
    fobj.close()

def toAzimuth(latevent,lonevent,latsource,lonsource):
    
        # Convert to radians.
        lat1 = math.radians(latsource);
        lon1 = math.radians(lonsource);
        lat2 = math.radians(latevent);
        lon2 = math.radians(lonevent);

       # Compute the angle.
        x = math.sin(lon1-lon2 )*math.cos(lat2);
        y = math.cos(lat1)*math.sin(lat2)-math.sin(lat1)*math.cos(lat2)*math.cos(lon1-lon2);

        angle = -math.atan2(x,y);

        if angle < 0.0 :
         angle  += math.pi * 2.0;

       #And convert result to degrees.
        angle = math.degrees(angle)
        angle = '%02f'%angle

        return angle;



#===============================================================================
def staltatriggering(sembmaxvalue,sembmaxlat,sembmaxlon,ntimes,Config,Folder):
    
    logger.info('%s' % ('Enter STA LTA TRIGGER') )
    mean2 = 0;
    ratio=1.5
    chkwin=100
    imax=0
    has_onset=0
    has_onset_dur=0

    fobjsembmax = open(os.path.join(Folder['semb'],'sembmaxvalue.txt'),'w')
    for i in range(ntimes):
        fobjsembmax.write('%d %d %f %.2f %.2f\n'%(i,i*int(Config['step']),sembmaxvalue[i],sembmaxlat[i],sembmaxlon[i]))
    fobjsembmax.close()    
# 
    fobjduration = open(os.path.join(Folder['semb'],'duration.txt'),'w')
    for k in range(ntimes):
        if k > 4 and k < 25:
            mean2 +=sembmaxvalue[k]

    mean2 /= 20.;
    lta = mean2;
    n = ntimes;
    
    for i in range(15, n-2):
        sta = (sembmaxvalue[i+1]+sembmaxvalue[i]+sembmaxvalue[i+2])/3
        if has_onset == 0:
            if sta/lta > ratio:
                ionset = 1
                fobjduration.write('# Initial onset at i = %d time = %d  semblance = %f sta: %f lta: %f sta/lta: %f\n'%(i,i*int(Config['step']),sembmaxvalue[i],sta,lta,sta/lta))
#                start = i*int(Config['step'])
                has_onset = 1
        else:
            if sembmaxvalue[i] > max and i <= i+chkwin:
#                max = sembmaxvalue[i]
                imax = i


    ratio = sembmaxvalue[imax]/lta;

    if ratio >= 1000.0:
        semb_start = sembmaxvalue[imax]/2.0
        semb_end   = sembmaxvalue[imax]/(ratio*0.01)
    elif ratio >= 100.0 and ratio < 1000.0:
        semb_start = sembmaxvalue[imax]/2.0
        semb_end   = sembmaxvalue[imax]/(ratio*0.1)
    elif ratio >= 5.0 and ratio < 100.0:
        semb_start = sembmaxvalue[imax]/2.5
        semb_end   = sembmaxvalue[imax]/1.5
    elif ratio >= 3.0 and ratio < 5.0:
        semb_start = sembmaxvalue[imax]/1.5
        semb_end   = sembmaxvalue[imax]/1.5
    else:
        semb_start = sembmaxvalue[imax]/1.3
        semb_end   = sembmaxvalue[imax]/(ratio*0.6)

#   search for duration as the time where semblance is above a portion of the semblance maximum
    for i in range (ionset,n-1):
        if has_onset_dur == 0: 
            if sembmaxvalue[i] >= semb_start:
                ionset_dur = i
                fobjduration.write('# Shifted event onset at i = %d time= %d semblance= %f\n'%(i,i*int(Config['step']),sembmaxvalue[i]))
                sonset = i*int(Config['step'])

                has_onset_dur = 1;
                fobjduration.write('%d %d %f %.2f %.2f\n' % (i,i*int(Config['step']),sembmaxvalue[i],sembmaxlat[i-ionset],sembmaxlon[i-ionset]))
        else:
            fobjduration.write('%d %d %f %.3f %.3f\n' % (i,i*int(Config['step']),sembmaxvalue[i],sembmaxlat[i-ionset],sembmaxlon[i-ionset]))
# 
#                                // duration:
            if sembmaxvalue[i] <= semb_end and i>=imax:
                iend_dur = i
                ndata = iend_dur - ionset_dur + 1
                fobjduration.write('# Shiftet event stop at i = %d time = %d semblance = %f\n'% (i,i*int(Config['step']),sembmaxvalue[i]))
                fobjduration.write('# Event duration = %d '%(i*int(Config['step'])-sonset))
                has_onset_dur=0
                break

    fobjduration.close()



def doCalculation(Config,WaveformDict,FilterMetaData,mint,TTTGridMap,Folder,Origin):
    
    logger.info('%s' % ('Enter Semblance Calculation') )
    
    new_frequence = int(Config['new_frequence'])
    winlen        = int(Config['winlen'])
    step          = float(Config['step'])
    forerun       = int(Config['forerun'])
    duration      = int(Config['duration'])
    dimX          = int(Config['dimX'])
    dimY          = int(Config['dimY'])
    gridspacing   = float(Config['gridspacing'])
    nostat = len(WaveformDict)
    traveltimes   = {}
    recordstarttime = ''
    minSampleCount = 999999999
    usedarrays = 5
    ntimes = int((forerun + duration)/step)
    nsamp = int(winlen * new_frequence);        #Anzahl der Samples pro Zeitschritt
    nstep = int(step * new_frequence);          #number of samples between subse
    
    
    ############################################################################
    calcStreamMap = WaveformDict
    for trace in calcStreamMap.iterkeys():
        recordstarttime = calcStreamMap[trace][0].stats.starttime
        d = calcStreamMap[trace][0].stats.starttime
        d = d.timestamp
        if calcStreamMap[trace][0].stats.npts < minSampleCount:
            minSampleCount = calcStreamMap[trace][0].stats.npts
    ############################################################################
    traces = np.ndarray(shape=(len(calcStreamMap),minSampleCount+1),dtype=float)
    traveltime = np.ndarray(shape=(len(calcStreamMap),dimX*dimY),dtype=float)
    latv = np.ndarray(dimX*dimY,dtype=float)
    lonv = np.ndarray(dimX*dimY,dtype=float)
    ############################################################################

    #for i in TTTGridMap.iterkeys():
     #   for j in TTTGridMap[i].GridArray.iterkeys():
      #      print j,TTTGridMap[i].GridArray[j].tt
     
    
    c=0
    streamCounter = 0
    print 'MSC: ',minSampleCount
    for key in calcStreamMap.iterkeys():
        
        streamID = key
        c2 = 0
        for i in calcStreamMap[key][0]:
            traces[c][c2] = i
            #print 'C: ',c,' C2: ',c2
            c2 += 1
        
        for key in TTTGridMap.iterkeys():
            
            if streamID == key:
                print "IN BEIDEN DA", streamID, key
                traveltimes[streamCounter] = TTTGridMap[key]
            else:
                "NEIN", streamID, key
        
        g = traveltimes[streamCounter]
        dimZ = g.dimZ
        mint = g.mint
        maxt = g.maxt
        Latul = g.Latul
        Lonul = g.Lonul
        Lator = g.Lator
        Lonor = g.Lonor
        
        gridElem = g.GridArray
        
        for x in range(dimX):
            for y in range(dimY):
                elem = gridElem[x, y]
                
                traveltime[c][x * dimY + y] = elem.tt
                latv[x * dimY + y] = elem.lat
                lonv[x * dimY + y] = elem.lon

        c += 1
        streamCounter += 1
    
        
    sembDict = {}
    sembmaxvalue = []
    sembmaxlat   = []
    sembmaxlon   = []
    
    rc = UTCDateTime(Origin['time'])
    rcs = '%s-%s-%s_%02d:%02d:%02d'% (rc.day,rc.month,rc.year, rc.hour,rc.minute,rc.second)
    
    
    fobjsembmax = open(os.path.join(Folder['semb'],'sembmax.txt'),'w')
    for i in range(ntimes):

        logger.info('Zeitschritt %d' % i)
        
        fobj = open(os.path.join(Folder['semb'],'%03d.ASC' % i),'w')
        fobj.write('# %s , %s\n' % (d,rcs))
        fobj.write('# step %ds| ntimes %d| winlen: %ds\n' % (step,ntimes,winlen))
        fobj.write('# \n')
        fobj.write('# southwestlat: %.2f dlat: %f nlat: %f \n'%(Latul,gridspacing,dimX))
        fobj.write('# southwestlon: %.2f dlon: %f nlat: %f \n'%(Lonul,gridspacing,dimY))
        fobj.write('# ddepth: 0 ndepth: 1 \n')
        
        sembmax = 0
        sembmaxX = 0
        sembmaxY = 0
        
        for j in range(dimX):
            
            for m in range(dimY):
                nomin=0
                denom=0
                semb = 0
                
                for l in range(nsamp):
                       summe = 0
                
                       for k in range(nostat):
                            relstart_samples = (int)((traveltime[k][j*dimY+m]-mint) * new_frequence + 0.5) + i*nstep;
                            #summe += traces[k][relstart_samples+l]
                            #denom += traces[k][relstart_samples+l]*trace[k][relstart_samples+l]
                            tmp = traces[k][relstart_samples + l]
                            summe += tmp        
                            denom += tmp ** 2
                    
                       nomin += summe*summe;
                    
                x = latv[j*dimY+m]
                y = lonv[j*dimY+m]
                semb = nomin/(float(nostat)*denom);

                fobj.write('%.2f %.2f %f\n' % (x,y,semb))
                
                if  semb > sembmax:
                       sembmax = semb;# search for maximum and position of maximum on semblance grid for given time step         
                       sembmaxX = latv[j*dimY+m];
                       sembmaxY = lonv[j*dimY+m];
        
        
        delta = locations2degrees(float(sembmaxX), float(sembmaxY), float(Origin['lat']), float(Origin['lon']))
        azi = toAzimuth(float(Origin['lat']), float(Origin['lon']),float(sembmaxX), float(sembmaxY))
        #print i,sembmax,sembmaxX,sembmaxY,' DELTA: ',delta,' OLAT: ',Origin['lat'],' OLON: ',Origin['lon'],float(sembmaxX), float(sembmaxY)
        sembmaxvalue.append(sembmax)
        sembmaxlat.append(sembmaxX)
        sembmaxlon.append(sembmaxY)
        
        fobjsembmax.write('%d %.2f %.2f %f %d %03f %f %03f\n' % (i*step,sembmaxX,sembmaxY,sembmax,usedarrays,delta,float(azi),delta*119.19))
        
        fobj.close()

    fobjsembmax.close()
    sembDict['value'] = sembmaxvalue
    sembDict['lat'] = sembmaxlat
    sembDict['lon'] = sembmaxlon
    #staltatriggering(sembmaxvalue,sembmaxlat,sembmaxlon,ntimes)
    
    return sembDict




def run(eventpath):

    Origin = readEvent(eventpath)
    Config = readConfig(eventpath)
    Meta   = readMetaInfoFile(eventpath)
    Folder = createFolder(eventpath)
    cpSkeleton(Folder,Config)
    writeConfig(Config,Origin,Folder)
    ntimes = int((int(Config['forerun']) + int(Config['duration']))/int(Config['step']))
    print 'NTIMES ',ntimes

    ne = Config['networks'].split(',')
    flag = 1
    for i in ne:
        network = Config[i].split(',')
        #print network

        FilterMeta = filterStations(Meta,Config,Origin,network)
        mint,maxt,TTTGridMap = calcTTT(Config,FilterMeta,Origin)
        tw = calculateTimeWindows(mint,maxt,Config,Origin)
        stream,Wd = readWaveforms(FilterMeta,tw,eventpath,Origin)
        Wdf = process(Wd,Config,Folder,i,FilterMeta)
        writeStationFile(FilterMeta,Folder,flag)
        t1 = time.time()
        sD = doCalculation(Config,Wdf,FilterMeta,mint,TTTGridMap,Folder,Origin)
        t2 = time.time()
        staltatriggering(sD['value'],sD['lat'],sD['lon'],ntimes,Config,Folder)
        flag +=1
        print '%s took %0.3f s' % ('Calculation', (t2-t1))


if __name__ == "__main__":
    print evpath
    run(evpath)