import os
import sys
import platform

import math
from   math     import radians, cos, sin, atan2
import logging
import time

# add local directories to import path

sys.path.append ('../tools/')                     
sys.path.append ('../Common/')
 
#
from  obspy.core.utcdatetime import UTCDateTime

#    Import from NumPy
import numpy as np          
import obspy
#    Import from common

import  Basic
import  Globals
import  Logfile
import  DataTypes
from    DataTypes  import Location
from    ObspyFkt   import loc2degrees 
from    ConfigFile import ConfigObj, OriginCfg 

import time
import os
import scipy
import numpy as num
import matplotlib.pyplot as plt
import plotly.plotly as py
from collections import OrderedDict



from pyrocko import gf, trace
from pyrocko import moment_tensor as mtm
from pyrocko.gf import ws, LocalEngine, Target, DCSource
from pyrocko import util, pile, model, config, trace, io, pile, catalog

km = 1000.


#    Import from Process

import  trigger


USE_C_CODE  =  False                   #16.12.2015

if USE_C_CODE : 
   import  CTrig                      
   import  Cm                          

else : 
   from semp import otest      

# -------------------------------------------------------------------------------------------------

logger = logging.getLogger('ARRAY-MP')

class SembMax (object):
    '''
    class to store sembmax object for each grid point
    '''
    def __init__(self,lat,lon,semb):

        self.lat  = lat
        self.lon  = lon
        self.semb = semb

# -------------------------------------------------------------------------------------------------

class FileSembMax(object):
    '''
    class to strore sembmax object for the sembmaxvalue file
    '''
    def __init__(self,istep,sembmaxX,sembmaxY,sembmax,usedarrays,delta,azi,deltakm):

        self.istep      = istep
        self.sembmaxX   = sembmaxX
        self.sembmaxY   = sembmaxY
        self.sembmax    = sembmax
        self.usedarrays = usedarrays
        self.delta      = delta
        self.azi        = azi
        self.deltakm    = deltakm

    def get(self):
        return ('%d %.2f %.2f %f %d %03f %f %03f\n' % (self.istep,self.sembmaxX,self.sembmaxY,
                                                       self.sembmax,self.usedarrays,self.delta,
                                                       self.azi,self.delta*119.19))

# -------------------------------------------------------------------------------------------------

def toAzimuth (latevent,lonevent,latsource,lonsource):
        '''
        method to calculate azimuth between two points
        '''
        # Convert to radians.
        lat1 = radians (latsource);
        lon1 = radians (lonsource);
        lat2 = radians (latevent);
        lon2 = radians (lonevent);

       # Compute the angle.
        x     =  sin(lon1-lon2 ) * cos(lat2);
        y     =  cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1-lon2);
        angle = -atan2 (x,y);

        if angle < 0.0 :
         angle  += math.pi * 2.0;

       #And convert result to degrees.
        angle = math.degrees (angle)
        angle = '%02f'%angle

        return angle;
    
# -------------------------------------------------------------------------------------------------

def writeSembMatricesSingleArray (SembList,Config,Origin,arrayfolder,ntimes,switch):
    '''
    method to write semblance matrizes from one processes to file for each timestep
    '''
    logger.info ('start write semblance matrices')
    
    cfg    = ConfigObj (dict=Config)
    origin = OriginCfg (Origin)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    latv   = []
    lonv   = []

    gridspacing = cfg.Float ('gridspacing')
    migpoints   = dimX * dimY

    o_lat       = origin.lat()         # float (Origin['lat'])
    o_lon       = origin.lon()         # float (Origin['lon'])
    oLatul      = 0
    oLonul      = 0
    
    z=0

    for i in xrange(dimX):
         oLatul = o_lat - ((dimX-1)/2) * gridspacing + i*gridspacing

         if z == 0 and i == 0:
             Latul = oLatul
         o=0    

         for j in xrange (dimY):
               oLonul = o_lon - ((dimY-1)/2) * gridspacing + j * gridspacing
               
               if o==0 and j==0:  Lonul = oLonul
               
               latv.append (oLatul)
               lonv.append (oLonul)
    #endfor

    rc  = UTCDateTime (Origin['time'])
    rcs = '%s-%s-%s_%02d:%02d:%02d'% (rc.day,rc.month,rc.year, rc.hour,rc.minute,rc.second)
    d   = rc.timestamp
    
    for a, i in enumerate(SembList):
        #logger.info('timestep %d' % a)

        fobj = open (os.path.join (arrayfolder,'%s-%s_%03d.ASC' % (switch,Origin['depth'],a)),'w')
        fobj.write ('# %s , %s\n' % (d,rcs))
        fobj.write ('# step %ds| ntimes %d| winlen: %ds\n' % (step,ntimes,winlen))
        fobj.write ('# \n')
        fobj.write ('# southwestlat: %.2f dlat: %f nlat: %f \n'%(Latul,gridspacing,dimX))
        fobj.write ('# southwestlon: %.2f dlon: %f nlon: %f \n'%(Lonul,gridspacing,dimY))
        fobj.write ('# ddepth: 0 ndepth: 1 \n')
        
        for j in range (migpoints):
            x    = latv[j]
            y    = lonv[j]
            z    = origin.depth()         # float(Origin['depth'])
            semb = i[j]

            fobj.write ('%.2f %.2f %.2f %.20f\n' % (x,y,z,semb))
        #endfor
            
        fobj.close()
    #endfor
         
# -------------------------------------------------------------------------------------------------

def collectSemb (SembList,Config,Origin,Folder,ntimes,arrays,switch):
    '''
    method to collect semblance matrizes from all processes and write them to file for each timestep
    '''
    Logfile.add ('start collect in collectSemb')
    
    cfg    = ConfigObj (dict=Config)
    origin = ConfigObj (dict=Origin)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    latv        = []
    lonv        = []

    gridspacing = cfg.Float ('gridspacing')
    migpoints   = dimX * dimY
    o_lat       = origin.lat()         # float (Origin['lat'])
    o_lon       = origin.lon()         # float (Origin['lon'])
    oLatul      = 0
    oLonul      = 0
    
    z=0

    for i in xrange(dimX):
         oLatul = o_lat - ((dimX-1)/2) * gridspacing + i*gridspacing

         if z == 0 and i == 0 :
             Latul = oLatul
         o=0    

         for j in xrange (dimY):
               oLonul = o_lon - ((dimY-1)/2) * gridspacing + j*gridspacing
               
               if o==0 and j==0:
                    Lonul = oLonul
               
               latv.append (oLatul)
               lonv.append (oLonul)
    #endfor  
  
    #print 'SL: ',SembList, type(SembList),type(SembList[0]),SembList[0].ndim

    tmp=1
    for a in SembList:
        tmp *= a

    #sys.exit()
    
    sembmaxvaluev = np.ndarray (ntimes,dtype=float)
    sembmaxlatv   = np.ndarray (ntimes,dtype=float)
    sembmaxlonv   = np.ndarray (ntimes,dtype=float)
    
    rc         = UTCDateTime(Origin['time'])
    rcs        = '%s-%s-%s_%02d:%02d:%02d'% (rc.day,rc.month,rc.year, rc.hour,rc.minute,rc.second)
    d          = rc.timestamp
    usedarrays = 5
    
    folder      = Folder['semb']
    fobjsembmax = open (os.path.join (folder,'sembmax_%s.txt' % (switch)),'w')
    
    for a, i in enumerate(tmp):
        logger.info('timestep %d' % a)
        

        fobj  = open (os.path.join (folder,'%s-%s_%03d.ASC' % (switch,Origin['depth'],a)),'w')
        #fobj = open (os.path.join (folder, '%03d.ASC'    % a),'w')

        fobj.write ('# %s , %s\n' % (d,rcs))
        fobj.write ('# step %ds| ntimes %d| winlen: %ds\n' % (step,ntimes,winlen))
        fobj.write ('# \n')
        fobj.write ('# southwestlat: %.2f dlat: %f nlat: %f \n'%(Latul,gridspacing,dimX))
        fobj.write ('# southwestlon: %.2f dlon: %f nlon: %f \n'%(Lonul,gridspacing,dimY))
        fobj.write ('# ddepth: 0 ndepth: 1 \n')
        
        
        sembmax  = 0
        sembmaxX = 0
        sembmaxY = 0
        
        origin = DataTypes.dictToLocation (Origin)
        uncert = np.std(i) #maybe not std?
        for j in range(migpoints):
            x    = latv[j]
            y    = lonv[j]
            semb = i[j]
  
            fobj.write ('%.2f %.2f %.20f\n' % (x,y,semb))
            
            if  semb > sembmax:
                sembmax  = semb;# search for maximum and position of maximum on semblance grid for given time step         
                sembmaxX = x;
                sembmaxY = y;

        delta = loc2degrees (Location (sembmaxX, sembmaxY), origin)
        azi   = toAzimuth   (float(Origin['lat']), float(Origin['lon']),float(sembmaxX), float(sembmaxY))

        sembmaxvaluev[a] = sembmax
        sembmaxlatv[a]   = sembmaxX
        sembmaxlonv[a]   = sembmaxY
        
        fobjsembmax.write ('%d %.2f %.2f %.20f %.20f %d %03f %f %03f\n' % (a*step,sembmaxX,sembmaxY,sembmax,uncert,usedarrays,delta,float(azi),delta*119.19))
        fobj.close()


    fobjsembmax.close()

    durationpath  = os.path.join (folder, "duration.txt")
    trigger.writeSembMaxValue (sembmaxvaluev,sembmaxlatv,sembmaxlonv,ntimes,Config,Folder)
    print 'DD2: ',durationpath
    trigger.semblancestalta (sembmaxvaluev,sembmaxlatv,sembmaxlonv)                          #hs
# -------------------------------------------------------------------------------------------------

def toMatrix (npVector, nColumns) :

    t   = npVector.tolist()[0]
    n   = nColumns
    mat = []

    for i in range (len(t) / n) :
       pos1  = i * n
       pos2  = pos1 + n
       slice = t [pos1:pos2]
       assert len(slice) == nColumns
       mat.append (slice)

    return mat




def  doCalc (flag,Config,WaveformDict,FilterMetaData,Gmint,Gmaxt,TTTGridMap,Folder,Origin, ntimes, switch, ev,arrayfolder):
    '''
    method for calculating semblance of one station array
    '''
    Logfile.add ('PROCESS %d %s' % (flag,' Enters Semblance Calculation') )
    Logfile.add ('MINT  : %f  MAXT: %f Traveltime' % (Gmint,Gmaxt))

    cfg = ConfigObj (dict=Config)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    new_frequence   = cfg.newFrequency()          #  ('new_frequence')
    forerun         = cfg.Int  ('forerun')
    duration        = cfg.Int  ('duration')
    gridspacing     = cfg.Float('gridspacing')

    nostat          = len (WaveformDict)
    traveltimes     = {}
    recordstarttime = ''
    minSampleCount  = 999999999

    ntimes = int ((forerun + duration)/step)
    nsamp  = int (winlen * new_frequence)
    nstep  = int (step   * new_frequence)
    from pyrocko import obspy_compat
    from pyrocko import trace as troll
    from pyrocko import orthodrome, model
    obspy_compat.plant()
    from pprint import pprint
    ############################################################################
    calcStreamMap = WaveformDict


    stations = []
    py_trs = []
    for trace in calcStreamMap.iterkeys():

        py_tr = obspy_compat.to_pyrocko_trace(calcStreamMap[trace])
        py_trs.append(py_tr)
        for il in FilterMetaData:
		if str(il) == str(trace):
                        szo = model.Station(lat=il.lat, lon=il.lon, station=il.sta, network=il.net, channels='Z', elevation=il.ele, location=il.loc)
			stations.append(szo)



    syn = True
    if syn == True:

	from pyrocko.gf import LocalEngine, Target, RectangularSource, ws
	from pyrocko import gf
	from pyrocko import trace as troll
        from pyrocko import pile
	from pyrocko.marker import PhaseMarker
	from beam_stack import BeamForming
	from pyrocko.gf import meta
	pjoin = os.path.join
	from abedeto.request import DataProvider, CakeTiming

	store_id = 'global_2hz'
	#engine = LocalEngine(store_superdirs=['/media/asteinbe/data/asteinbe/aragorn/andreas/Tibet/'])
	engine = LocalEngine(store_superdirs=['/media/asteinbe/memory47/'])


	bounds = OrderedDict([
	    ('north_shift', (-20.*km, 20.*km)),
	    ('east_shift', (-20.*km, 20.*km)),
	    ('depth', (3.*km, 8.*km)),
	    ('magnitude', (6.2, 6.4)),
	    ('strike', (100., 140.)),
	    ('dip', (40., 60.)),
	    ('rake', (-100, -150.)),
	    ('timeshift', (-20., 20.)),
	    ])



        targets =[]
	for st in stations:
		target= Target(
		lat=st.lat,
		lon=st.lon,
		store_id=store_id,
		codes=(st.network, st.station, st.location, 'BHZ'))
		targets.append(target)


        ## scaling law?
	timeev = util.str_to_time('2016-11-25 14:24:30.000')
	source_dc = RectangularSource(
	    lat=float(ev.lat),
	    lon=float(ev.lon),
	    depth=ev.depth*1000.,
	    strike=110.,
	    dip=80.,
	    rake=-160.,
	    width=10.*1000.,
	    length=55.*1000.,
	    nucleation_x=0.,
       slip = 0.8,
            nucleation_y=0.0,
	    #stf=gf.BoxcarSTF(duration=20.0),
	    time = timeev,)
 	    #magnitude=6.7)

	response = engine.process(source_dc, targets)
           
	synthetic_traces = response.pyrocko_traces()

	i =0
	trs_org= []
	trs_orgs= []
	fobj = os.path.join (arrayfolder,'shift.dat')
	xy = num.loadtxt(fobj, usecols=1, delimiter=',')
   #import csv
	#with open(fobj, "w") as output:
    #    	writer = csv.writer(fobj, lineterminator='\n')
     #   	writer.writerows(res)

        calcStreamMapsyn= calcStreamMap.copy()
        for trace in calcStreamMapsyn.iterkeys():
		mod = synthetic_traces[i]
        	recordstarttime = calcStreamMapsyn[trace].stats.starttime.timestamp
        	recordendtime = calcStreamMapsyn[trace].stats.endtime.timestamp
		#tmin_new = util.str_to_time('2009-04-06 01:32:42.000')
		mod.shift(recordstarttime-mod.tmin)
		extracted = mod.chop(recordstarttime, recordendtime, inplace=False)
		tr_org = obspy_compat.to_pyrocko_trace(calcStreamMapsyn[trace])

	#	tr_org.shift(xy[i])


	#	tr_org.bandpass(4,0.025,0.24)
	#	tr_org.downsample_to(0.5)

    #        	displacement = tr_org.transfer(
     #           10.,  # rise and fall of time domain taper in [s]
      #          (0.005, 0.01, 1., 2.))
	#	trs_org.append(displacement)
		#troll.snuffle(tr_org)
		#calcStreamMap[trace]=obspy_compat.to_obspy_trace(tr_org)

		#mod.ydata = mod.ydata*0.
		#extracted.bandpass(4,0.025,0.24)
		#extracted.downsample_to(0.5)
        	synthetic_obs_tr = obspy_compat.to_obspy_trace(extracted)

		calcStreamMapsyn[trace]=synthetic_obs_tr
		trs_orgs.append(tr_org)

		trs_org.append(extracted)
		#troll.snuffle(extracted)
		i = i+1
    	#	calcStreamMap = calcStreamMapsyn
	#store = engine.get_store(store_id)
   
#	timing = meta.Timing('first(P)')
	timing = CakeTiming(
       phase_selection='first(p|P|PP|P(cmb)P(icb)P(icb)p(cmb)p)-20',
       fallback_time=100.)
	traces = trs_orgs
	event = model.Event(lat=float(ev.lat), lon=float(ev.lon), depth=ev.depth*1000., time=timeev)
	directory = arrayfolder
	bf = BeamForming(stations, traces, normalize=True)
	shifted_traces = bf.process(event=event,
              timing=timing,
              fn_dump_center=pjoin(directory, 'array_center.pf'),
              fn_beam=pjoin(directory, 'beam.mseed'))
	#bf.snuffle()
#	troll.snuffle(trs_orgs)
#	troll.snuffle(trs_org)
	calcStreamMapshifted= calcStreamMap.copy()
	i = 0
	for trace in calcStreamMapsyn.iterkeys():
         mod = shifted_traces[i]
         extracted = mod.chop(recordstarttime, recordendtime, inplace=False)
         shifted_obs_tr = obspy_compat.to_obspy_trace(extracted)
         calcStreamMapsyn[trace]=shifted_obs_tr    
         i = i+1
         
	calcStreamMap = calcStreamMapshifted

    for trace in calcStreamMap.iterkeys():
        recordstarttime = calcStreamMap[trace].stats.starttime
        d = calcStreamMap[trace].stats.starttime
        d = d.timestamp

        if calcStreamMap[trace].stats.npts < minSampleCount:
            minSampleCount = calcStreamMap[trace].stats.npts

    ############################################################################
    traces     = np.ndarray (shape=(len(calcStreamMap), minSampleCount), dtype=float)
    traveltime = np.ndarray (shape=(len(calcStreamMap), dimX*dimY), dtype=float)
    latv       = np.ndarray (dimX*dimY, dtype=float)
    lonv       = np.ndarray (dimX*dimY, dtype=float)
    ############################################################################


    c=0
    streamCounter = 0
    
    for key in calcStreamMap.iterkeys():
        streamID = key
        c2       = 0

        for o in calcStreamMap[key]:
            if c2 < minSampleCount:
                traces[c][c2] = o

                c2 += 1


        for key in TTTGridMap.iterkeys():
            
            if streamID == key:
                traveltimes[streamCounter] = TTTGridMap[key]
            else:
                "NEIN", streamID, key


        if not streamCounter in traveltimes : 
           continue                              #hs : thread crashed before

        g     = traveltimes[streamCounter]
        dimZ  = g.dimZ
        mint  = g.mint
        maxt  = g.maxt
        Latul = g.Latul
        Lonul = g.Lonul
        Lator = g.Lator
        Lonor = g.Lonor
        
        gridElem = g.GridArray
        
        for x in range(dimX):
            for y in range(dimY):
                elem = gridElem[x, y]
                
                traveltime [c][x * dimY + y] = elem.tt
                latv [x * dimY + y] = elem.lat
                lonv [x * dimY + y] = elem.lon
        #endfor

        c += 1
        streamCounter += 1
    #endfor


    ############################## CALCULATE PARAMETER FOR SEMBLANCE CALCULATION ##################
    nsamp     = winlen * new_frequence
    nstep     = int (step*new_frequence)
    migpoints = dimX * dimY

    dimZ = 0
    new_frequence = cfg.newFrequency ()              # ['new_frequence']
    maxp = int (Config['ncore'])


    Logfile.add ('PROCESS %d  NTIMES: %d' % (flag,ntimes))
  
    if False :
       print ('nostat ',nostat,type(nostat))
       print ('nsamp ',nsamp,type(nsamp))
       print ('ntimes ',ntimes,type(ntimes))
       print ('nstep ',nstep,type(nstep))
       print ('dimX ',dimX,type(dimX))
       print ('dimY ',dimY,type(dimY))
       print ('mint ',Gmint,type(mint))
       print ('new_freq ',new_frequence,type(new_frequence))
       print ('minSampleCount ',minSampleCount,type(minSampleCount))
       print ('latv ',latv,type(latv))
       print ('traces',traces,type(traces))
       print ('traveltime',traveltime,type(traveltime))


    try:
    	cs = cfg.cs ()
    except:
	cs = 1
    if cs == 1:
    	    csmaxvaluev = np.ndarray (ntimes,dtype=float)
    	    csmaxlatv   = np.ndarray (ntimes,dtype=float)
    	    csmaxlonv   = np.ndarray (ntimes,dtype=float)
	    folder      = Folder['semb']
            fobjcsmax = open (os.path.join (folder,'csmax_%s.txt' % (switch)),'w') 
	    traveltimes = traveltime.reshape (1,nostat*dimX*dimY)
	    traveltime2 = toMatrix (traveltimes, dimX * dimY)  # for relstart
	    import matplotlib as mpl
	    traveltime = traveltime.reshape (dimX*dimY,nostat)
	    import scipy.optimize as spopt
	    import scipy.fftpack as spfft
	    import scipy.ndimage as spimg
	    import cvxpy as cvx
	    import matplotlib.pyplot as plt
	    A = spfft.idct(traveltime, norm='ortho', axis=0)
	    #A = traveltime
	    n = (nostat*dimX*dimY)
	    vx = cvx.Variable(dimX*dimY)
	    res = cvx.Variable(1)
	    objective = cvx.Minimize(cvx.norm(res, 1))
	    back2 = np.zeros([dimX,dimY])
	    l = int(nsamp)
            fobj  = open (os.path.join (folder,'%s-%s_%03d.cs' % (switch,Origin['depth'],l)),'w')
        #fobj = open (os.path.join (folder, '%03d.ASC'    % a),'w')

            
	    for i in range (ntimes) :
	    		    ydata = []
                            try:
			    	    for tr in traces:
						relstart = int ((dimX*dimY - mint) * new_frequence + 0.5) + i * nstep
						tr=spfft.idct(tr[relstart+i:relstart+i+dimX*dimY], norm='ortho', axis=0)

						ydata.append(tr)
				    ydata = np.asarray(ydata)
				    ydata     = ydata.reshape     (dimX*dimY,nostat)
				    #print np.shape(A), np.shape(vx), np.shape(ydata)
				    
				    constraints = [res == cvx.sum_entries( 0+ np.sum([ydata[:,x]-A[:,x]*vx  for x in range(nostat) ]) ) ]
				   # constraints = [0 <= ydata[0,:]- A*vx]
				   # constraints = [A[0,:]*vx == ydata[0,:]]
				    prob = cvx.Problem(objective, constraints)
				    result = prob.solve(verbose=False, max_iters=200) 

				    x = np.array(vx.value)
				    x = np.squeeze(x)

				    back1    = x.reshape     (dimX,dimY)
				    sig = spfft.idct(x, norm='ortho', axis=0) 
				    back2 = back2 + back1
				    xs = np.array(res.value)
				    xs = np.squeeze(xs)
				    max_cs = np.max(back1)
				    idx = np.where(back1==back1.max())
				    csmaxvaluev[i] = max_cs
				    csmaxlatv[i]   = latv[idx[0]]
				    csmaxlonv[i]   = lonv[idx[1]]
				    fobj.write ('%.5f %.5f %.20f\n' % (latv[idx[0]],lonv[idx[1]],max_cs))
				    fobjcsmax.write ('%.5f %.5f %.20f\n' % (latv[idx[0]],lonv[idx[1]],max_cs))
                            except:
			            pass
	    fobj.close()
	    fobjcsmax.close()
    t1 = time.time()
    traces     = traces.reshape     (1,nostat*minSampleCount)
    traveltime = traveltime.reshape (1,nostat*dimX*dimY)

    if USE_C_CODE :
       k  = Cm.otest (maxp,nostat,nsamp,ntimes,nstep,dimX,dimY,Gmint,new_frequence,
                      minSampleCount,latv,lonv,traveltime,traces)
    else :
       k = otest (maxp,nostat,nsamp,ntimes,nstep,dimX,dimY,Gmint,new_frequence,
                  minSampleCount,latv,lonv,traveltime,traces)                       #hs

    t2 = time.time()
    
    Logfile.add ('%s took %0.3f s' % ('CALC:', (t2-t1)))

 
    partSemb = k

    partSemb  = partSemb.reshape (ntimes,migpoints)

    
    return partSemb

