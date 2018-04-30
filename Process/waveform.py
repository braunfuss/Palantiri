
import os
import sys
import logging
import fnmatch
import platform

WINDOWS = (platform.system() == 'Windows')

# add local directories to import path 

sys.path.append ('../tools/')                     
sys.path.append ('../Common/')

from obspy.core             import read
from obspy.core.stream      import Stream, Trace
from obspy.core.utcdatetime import UTCDateTime

import numpy

#       Import from common

import  Basic
import  Globals
import  Logfile
from    ConfigFile import ConfigObj, FilterCfg 

from    config import Station,Event     #       Import from Tools
import  ttt                             #       Import from Process

# -------------------------------------------------------------------------------------------------

def getStation (streamID,MetaDict):

    for i in MetaDict:
        if fnmatch.fnmatch (streamID,i.getName()) :  stobj = i

    return stobj

# -------------------------------------------------------------------------------------------------

def getGain(streamID,MetaDict):
    
    gain=1
    
    for i in MetaDict:
        if fnmatch.fnmatch(streamID,i.getName()):
            gain = i.gain
            
    return gain
# -------------------------------------------------------------------------------------------------

def readWaveforms (stationList,tw,EventPath,Origin):
    
    t2      = UTCDateTime (Origin.time)
    sdspath = os.path.join (EventPath,'data',str(t2.year))

    Wdict = {}
    
    for i in stationList:
        streamData = i.getName() + '.D.' + str(t2.year) + '.' + str("%03d" % t2.julday)
        entry      = os.path.join(sdspath,i.net,i.sta,i.comp+'.D',streamData)
        tdiff      = tw['end']-tw['start']
        
        try:
            st = read (entry,format="MSEED",starttime=tw['start'],endtime=tw['end'],nearest_sample=True)

        except:
            Logfile.error ('readWaveforms : File not found', entry)
            pass

        if len (st.get_gaps()) > 0:
            st.merge (method=0, fill_value='interpolate', interpolation_samples=0)
            
        print '--------------------------------------------------------------------------------'
        
        if len(st) > 0:
            trdiff = st[0].stats.endtime-st[0].stats.starttime
            
            totaldiff = abs(trdiff - tdiff)
            #print streamData,' TDIFF ',tdiff,' TRDIFF ',trdiff,' TOTALDIFF ',totaldiff
            
            if totaldiff < 1:
                Wdict[i.getName()] = st
               #stream += st
                Logfile.add (i.getName() + ' added to StreamList ')
        else:
           print ' OUT ',streamData
        
        print st
    #endfor        

    Logfile.red ('%d Streams added with available Data' % len(Wdict))
    return Wdict

# -------------------------------------------------------------------------------------------------

def writeWaveform (Folder, station, Stream, flag, network):
    
    if   flag == 'U':  s1 = 'unfiltered'
    elif flag == 'F':  s1 = 'filtered'
    elif flag == 'R':  s1 = 'resampled'

    else : Logfile.abort ('writeWaveform : Illegal flag <' + flag +'>')

    fname = ('%s-%s-%s.mseed') % (network,station,flag)
    name  = os.path.join (Folder['mseed'],fname)
    Stream.write (name, format='MSEED')

    Logfile.add  ('%s stream for station %s written '% (s1, station))
    
# -------------------------------------------------------------------------------------------------

def resampleWaveform_2 (Waveform, end_frequence):
    return resampleWaveform (Waveform, end_frequence)

def resampleWaveform (Waveform, end_frequence):
        
        Logfile.add ('enter resampling')
        
        new_frequence = end_frequence
        sampling_rate = Waveform.stats.sampling_rate
        if sampling_rate == end_frequence:
            return Waveform
        
#       if int (sampling_rate) in [20, 40, 80, 100] :            #hs : 16.10.2014
        if int (sampling_rate) in [20, 40, 50,80, 100, 125, 200] :
            oldsr  = sampling_rate
            factor = int (sampling_rate / new_frequence + 0.5)

            Logfile.add ('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))
            #Waveform.decimate (factor)
            Waveform.resample(new_frequence, window='hanning', no_filter=True)
            return Waveform
    
        else :
          if sampling_rate == 25:               #upsampling with factor = 4
             n1 = 2
             n2 = 4

  #        elif sampling_rate == 50:             #upsampling with factor = 2
   #          n1 = 1
    #         n2 = 2
          #endif

#          else : Logfile.abort ('Illegal sampling rate ' + str(sampling_rate))
 #         print "before data"
  #        data   = numpy.zeros    ((Waveform.stats.npts * n1))
   #       test   = numpy.fft.rfft (Waveform.data)
    #      zero   = numpy.fft.rfft (data)
     #     t      = numpy.concatenate((test,zero))
      #    tbdata = numpy.fft.irfft(t)
       #     
 #         tm = Trace(tbdata)
#
     #     tm.stats.network       = Waveform.stats.network
      #    tm.stats.station       = Waveform.stats.station
       #   tm.stats.sampling_rate = Waveform.stats.sampling_rate * n2
        #  tm.stats.channel       = Waveform.stats.channel
         # tm.stats.starttime     = UTCDateTime(Waveform.stats.starttime)
          #tm.stats._format       = 'MSEED'
        
          #oldsr  = tm.stats.sampling_rate
          #factor = int(tm.stats.sampling_rate / new_frequence + 0.5)
         # print "before deci"
          #Logfile.add ('downsampling %s from %d Hz to %d Hz with factor after upsampling %d' % (tm, oldsr, new_frequence, factor))
          Waveform.resample(new_frequence, window='hanning', no_filter=True)
          print "decimated"
          return tm
            
# -------------------------------------------------------------------------------------------------
#filtering

def filterWaveform_2 (Config, wafeform, station) : 
    return filterWaveform (Config, wafeform, station)     

def filterWaveform (Config, wafeform, station) :      

    cfg    = FilterCfg (self.Config)
    switch = cfg.filterswitch ()              # 'filterswitch'
            
    if switch == 1:
       which = 'bandpass'
       waveform.filter (which, freqmin   = cfg.flo(),           # ['flo']
                               freqmax   = cfg.fhi(),           # ['fhi']
                               corners   = cfg.ns(),            # ['ns']
                               zerophase = bool (Config['zph']))
    elif switch == 2:
      which = "lowpass"
      waveform.filter (which, freq      = cfg.l_fc(),             # ['l_fc']
                              corners   = cfg.l_ns(),             # ['l_ns']
                              zerophase = bool (Config['l_zph']))
    elif switch == 3:
      which = 'highpass'
      waveform.filter (which, freq     = cfg.h_fc(),             # ['h_fc']
                              corners  = cfg.h_ns(),             # ['h_ns']
                              zerophase= bool (Config['h_zph']))
    else : return None

    Logfile.add ('%s filtered stream for station %s '% (which, i))
    return wafeform

# -------------------------------------------------------------------------------------------------

def processWaveforms (WaveformDict,Config,Folder,network,MetaDict,Event,switch,Xcorr):
    
    Logfile.red ('Start Processing')
    
    cfg           = FilterCfg (Config)
    new_frequence = cfg.newFrequency()                #  ['new_frequence']

    vp,vs,rho = ttt.calcak135parameter (Event)
    
    for index,i in enumerate (WaveformDict):
            Logfile.add ('%s:%s -------------------------------------------------------------' % (index,i))
            
            if Config ['export_unfiltered'].capitalize() == 'True':
               writeWaveform (Folder,i,WaveformDict[i],'U',network)

            # TODO implement gain multiplication
            station = getStation (i,MetaDict)

            psign = 1
            
#           print 'MetaDict ',MetaDict,station,station.takeoff
            #needs to be implemented correctly

            if cfg.Int ('fm') == 1:
                print 'Event ',Event, Event.strike,Event.dip,Event.rake
                azi   = ttt.bearing (Event, station)
                bazi  = ttt.backazi (station, Event)
                psign = ttt.dubcup  (rho, vp,vs,Event.strike,Event.dip,Event.rake, azi, station.takeoff)

               #print 'Takeoff ',station.takeoff,' Azi ',azi,' Bazi ',bazi,' vp ',vp,' vs ',vs,' rho ',rho,' Psign ',psign
                msg = 'Takeoff ' + str(station.takeoff) + ' Azi ' + str(azi) + ' Bazi ' + str(bazi)
                msg += (' vp ' + str(vp) + ' vs ' + str(vs) + ' rho ' + str(rho) + ' Psign ' + str(psign))
                Logfile.add (msg)
            
            '''
            if int(Config['xcorr']) == 1:
                #print Xcorr[i].value
                if Xcorr[i].value > 0:   psign = 1
                else:                    psign = -1
            '''
            try:
            	Logfile.add (station.getName() + ' ' + station.lat + ' ' + station.lon + ' ' + station.gain + ' PSIGN: ' + str(psign))
            except:
		try: 
			print psign, station.getName(), station.lat, station.lon, station.gain
		except:
			pass
            if psign == -1:
                Logfile.add ('correcting polarisation for station %s ' % (i))
                -1.*WaveformDict[i][0].data
            
            #remove offset
         #   WaveformDict[i].detrend (type='demean')
            
            #gain correctur
            gain = float (station.gain)

            if gain == 0.0 or gain == -1.0 : gain = 1    #hs : gain missing in metafile
            
            WaveformDict[i][0].data * (1.0 / gain)
            
            #filtering
            #switch = cfg.filterswitch ()              # 'filterswitch'
            
            if switch == 0:
                Logfile.add ('bandpass filtered stream for station %s '% (i))

                WaveformDict[i].filter ('bandpass',
                                        freqmin   = cfg.flo(),           # ['flo']
                                        freqmax   = cfg.fhi(),           # ['fhi']
                                        corners   = cfg.ns(),            # ['ns']
                                        zerophase = bool (Config['zph']))
                
            elif switch == 1:
                Logfile.add ('bandpass filtered stream for station %s '% (i))

                WaveformDict[i].filter ('bandpass',
                                        freqmin   = cfg.flo2(),           # ['flo']
                                        freqmax   = cfg.fhi2(),           # ['fhi']
                                        corners   = cfg.ns2(),            # ['ns']
                                        zerophase = bool (Config['zph2']))

    #        elif switch == 2:
     #          Logfile.add ('lowpass filtered stream for station %s '% (i))

      #         WaveformDict[i].filter ("lowpass",
        #                              freq      = cfg.l_fc(),             # ['l_fc']
       #                               corners   = cfg.l_ns(),             # ['l_ns']
         #                             zerophase = bool (Config['l_zph']))
          #  elif switch == 3:
           #     Logfile.add ('highpass filtered stream for station %s '% (i))

            #    WaveformDict[i].filter ("highpass", 
             #                          freq     = cfg.h_fc(),             # ['h_fc']
              #                         corners  = cfg.h_ns(),             # ['h_ns']
               #                        zerophase= bool (Config['h_zph']))
            else :
               Logfile.add ('no filter set for station %s '% (i))

            if Config ['export_filtered'].capitalize() == 'True':
               writeWaveform (Folder,i,WaveformDict[i],'F',network)

            #resampling
            j = resampleWaveform (WaveformDict[i][0],new_frequence)
            WaveformDict[i] = j
            #Logfile.add (WaveformDict[i])

            if Config ['export_resampled'].capitalize() == 'True':
               writeWaveform (Folder,i,WaveformDict[i],'R',network)
    #endfor            
    
    return WaveformDict
