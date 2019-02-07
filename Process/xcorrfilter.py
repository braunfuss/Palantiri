
import os
import sys

# add local directories to import path      

sys.path.append ('../tools/')                     
sys.path.append ('../Common/')
       
import logging
import fnmatch
from   ConfigParser import SafeConfigParser

import obspy.core
from   obspy.core             import read
from   obspy.core.utcdatetime import UTCDateTime
import obspy.signal.cross_correlation

#    Import from common

import  Basic
import  Globals
import  Logfile
from    ObspyFkt   import loc2degrees, obs_TravelTimes

logger = logging.getLogger('ARRAY-MP')

# --------------------------------------------------------------------------------------------------
class Xcorr(object):
    
    def __init__(self,Origin,StationMeta,EventPath):

        self.Origin  = Origin
        self.StationMeta = StationMeta
        self.EventPath   = EventPath

    # ----------------------------------------------------------------------------------------------
    
    def traveltimes(self):
        
        Logfile.red ('Enter AUTOMATIC FILTER')
        T = []
        Wdict = {}
        
        for i in self.StationMeta:        
            de = loc2degrees   (self.Origin, i)
            tt = obs_TravelTimes (de, self.Origin.depth)

            if tt[0]['phase_name'] == 'P':
               time = tt[0]['time']
               T.append(time)
        
            #print i,i.getName(),i.lat,i.lon,time
        
            tw = self.calculateTimeWindows (time)
            w  = self.readWaveformsCross (i,tw)
            Wdict [i.getName()] = w
        #endfor

        Logfile.red ('Exit AUTOMATIC FILTER')

        return Wdict
    # ----------------------------------------------------------------------------------------------
    
    def calculateTimeWindows(self,mint):
        tw = {}
        st = str(self.Origin.time)[:-1]

        tw['start'] = UTCDateTime(UTCDateTime(st)+(mint-5))
        tw['end']   = tw['start'] + 20

        Logfile.add (' ORIGIN TIME %s' % UTCDateTime (self.Origin.time))
        Logfile.add (' TIME WINDOW: %s - %s' % (tw['start'], tw['end']))
        return tw
    # ----------------------------------------------------------------------------------------------
    
    def filterWaveform(self,Waveform):
        
        Logfile.red ('Filter Waveform:')

        for i in Waveform:
            i.detrend ("simple")
            i.filter("bandpass", freqmin=0.05, freqmax=1, corners=3, zerophase=False)

        return Waveform
    # ----------------------------------------------------------------------------------------------
    
    def readWaveformsCross (self,station,tw):
        
        time = self.Origin.time
        ts   = time.split('T')
    
        datet  = ts[0]
        datet  = datet.split('-')
        year   = datet[0].strip()
        month  = datet[1]
        day= datet[2]
        #timep = ts[1][:-1]
     
        julday  = UTCDateTime (int(year),int(month),int(day)).julday
        julday  = "%03d" % julday
        sdspath = os.path.join (self.EventPath,'data',year)

        if station.loc =='--':
            station.loc =''

        streamData = station.net+'.'+station.sta+'.'+station.loc+'.'+station.comp+'.D.'+str(year)+'.'+str(julday)
        entry  = os.path.join (sdspath,station.net,station.sta,station.comp+'.D',streamData)
        #print entry
        st = read (entry,format="MSEED", starttime=tw['start'], endtime=tw['end'], nearest_sample=True)
        #print st
        
        if len (st.getGaps()) > 0:
            st.merge (method=0, fill_value='interpolate', interpolation_samples=0)

        st[0].stats.starttime = UTCDateTime(1000)

        stream = self.filterWaveform(st)
        return stream
    # ----------------------------------------------------------------------------------------------

    def doXcorr (self):

        StreamDict = self.traveltimes()
        corrDict   = {}

        for stream in StreamDict.iterkeys():
            ref = StreamDict[stream][0].data
    
        Logfile.red ('Enter Xcorr Procedure')

        for stream in StreamDict.iterkeys():
            #print stream, StreamDict[stream][0]
            a,b = obspy.signal.cross_correlation.xcorr (ref, StreamDict[stream][0], 100)
            shift   = a / StreamDict[stream][0].stats.sampling_rate
            corrDict[stream] = shift
            #print 'Index: ',a, ' Value: ',b,' ----> ',stream

        Logfile.red ('Leave Xcorr Procedure')
        return corrDict
