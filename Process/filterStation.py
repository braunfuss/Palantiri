
import os
import sys
import logging
import platform

WINDOWS =(platform.system() == 'Windows')

# add local directories to import path

sys.path.append('../tools/')
sys.path.append('../Common/')

from  obspy.core.stream       import Stream
from  obspy.core              import read
from  obspy.core.utcdatetime  import UTCDateTime

#       Import from Common

import  Basic
import  Globals
import  Logfile
import  Debug
from    ObspyFkt   import loc2degrees, obs_TravelTimes
from    ConfigFile import ConfigObj

#       Import from Tools

from config import Station

#       Import from Process

import  waveform
import  times

# --------------------------------------------------------------------------------------------------

def readWaveformsCross(station, tw, EventPath, Origin):

    time = Origin.time
    ts   = time.split('T')

    datet = ts[0]
    datet = datet.split('-')
    year  = datet[0].strip()
    month = datet[1]
    day   = datet[2]
   #timep = ts[1][:-1]

    julday  = UTCDateTime(int(year),int(month),int(day)).julday
    julday  = "%03d" % julday
    sdspath = os.path.join(EventPath,'data',year)

    #Wdict = {}

    streamData = station.getName()+'.D.'+str(year)+'.'+str(julday)
    entry      = os.path.join(sdspath,station.net,station.sta,station.comp+'.D',streamData)
    st         = read(entry,format="MSEED", starttime=tw['start'], endtime=tw['end'], nearest_sample=True)

    if len(st.getGaps()) > 0:
        st.merge(method=0, fill_value='interpolate', interpolation_samples=0)

    #Wdict[i.getName()] = st
    stream = st

    return stream
# --------------------------------------------------------------------------------------------------

def traveltimes(MetaDict,Config,Event,Folder,evpath):

    Logfile.red('Enter AUTOMATIC FILTER')
    T = []

    for i in MetaDict:
        delta = loc2degrees    (Event, i)
        tt    = obs_TravelTimes(delta, Event.depth)

        if tt[0]['phase_name'] == Config['ttphase'] :
           time = tt[0]['time']
           T.append(time)

        mint = min(T)
        maxt = max(T)
        
        tw = times.calculateTimeWindows(mint,maxt,Config,Event)
        readWaveformsCross(i,tw,evpath,Event)
    #endfor

    Logfile.red('Exit AUTOMATIC FILTER')
