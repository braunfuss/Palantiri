import config
import sys
from optparse import OptionParser
import logging
from obspy.signal import rotate
import fnmatch
import os
from obspy.core import read
from obspy.core import util
from obspy.core.stream import Trace
from obspy.core.utcdatetime import UTCDateTime

logger = logging.getLogger(sys.argv[0])

parser = OptionParser(usage="%prog -p Eventpath")
parser.add_option("-p", type="string", dest="evpath", help="evpath")

(options, args) = parser.parse_args()
evpath = options.evpath

C = config.Config(evpath)
Origin = C.parseConfig('origin')
cfg = ConfigObj (dict=C)
if cfg.pyrocko_download == True:
    Meta = C.readpyrockostations()
else:
    Meta = C.readMetaInfoFile()

print Origin
o_lat = float(Origin['lat'])
o_lon = float(Origin['lon'])

class StatRot(object):
    def __init__(self,sname,bazi,cl,ncomp='',ecomp=''):
        self.sname = sname
        self.bazi = bazi
        self.cl = cl
        self.ncomp = ncomp
        self.ecomp = ecomp

def create_dir(wdir):
  if os.access(wdir, os.W_OK):
    return True

  try:
    os.makedirs(wdir)
    return True
  except:
    return False


toDoList = []
sDict = {}

P = []
for root,dirs,files in os.walk(os.path.join(evpath,'data')):
        L = []
        for i in files:
            t = os.path.join(root,i)
            r = os.path.join('/',*t.split('/')[:-1])
            if not fnmatch.fnmatch(os.path.basename(r), '*HZ.D'):
                P.append(t)


T = []
for index,station in enumerate(Meta):
    if fnmatch.fnmatch(station.comp, '*HZ'):
        print 'INDEX ',index,station,station.lat,station.lon
        s_lat = float(station.lat)
        s_lon = float(station.lon)
        dist,azi,bazi = util.geodetics.gps2DistAzimuth(s_lat,s_lon,o_lat,o_lon)
        shortname = ('%s.%s')%(station.net,station.sta)
        for j in P:
            sdsname = os.path.basename(j)[:-11]
            ss = sdsname.split('.')
            sdsshort = ('%s.%s')%(ss[0],ss[1])
            if shortname == sdsshort:
                print 'JUHU',station, j
                st = read(j)
                tr = st[0]
                if fnmatch.fnmatch(tr.stats.channel,'*E'):
                    e = tr.data
                    print 'E COMP'
                if fnmatch.fnmatch(tr.stats.channel,'*N'):
                    n = tr.data
                    print 'N COMP'

        if len(n) == len(e):

            z = rotate.rotate_NE_RT(n, e, bazi)#return r and t data array


            tmr = Trace(z[0])
            tmr.stats.network = tr.stats.network
            tmr.stats.station = tr.stats.station
            tmr.stats.sampling_rate = tr.stats.sampling_rate
            tmr.stats.channel = 'R'
            tmr.stats.starttime = UTCDateTime(tr.stats.starttime)
            tmr.stats._format = 'MSEED'

            pr = os.path.join(os.path.join('/',*j.split('/')[:-2]),'R.D')
            os.makedirs(pr)
            print 'PR ',pr
            jd = "%03d" % tmr.stats.starttime.julday
            n = ('%s.%s.%s.%s.%s.%s')%(tmr.stats.network,tmr.stats.station,tmr.stats.location,tmr.stats.channel,str(tmr.stats.starttime.year),str(jd))
            tmr.write(os.path.join(pr,n),format='MSEED')

            tmt = Trace(z[1])
            tmt.stats.network = tr.stats.network
            tmt.stats.station = tr.stats.station
            tmt.stats.sampling_rate = tr.stats.sampling_rate
            tmt.stats.channel = 'T'
            tmt.stats.starttime = UTCDateTime(tr.stats.starttime)
            tmt.stats._format = 'MSEED'

            pt = os.path.join(os.path.join('/',*j.split('/')[:-2]),'T.D')
            os.makedirs(pt)
            print 'PT ',pt
            jd = "%03d" % tmr.stats.starttime.julday
            n = ('%s.%s.%s.%s.%s.%s')%(tmr.stats.network,tmr.stats.station,tmr.stats.location,tmr.stats.channel,str(tmr.stats.starttime.year),str(jd))
            tmt.write(os.path.join(pt,n),format='MSEED')


            print z

        print '\n\n\n\n'




'''
for i in T:
    L = []
    for j in P:
        sdsname = os.path.basename(j)[:-11]
        ss = sdsname.split('.')
        sdsshort = ('%s.%s')%(ss[0],ss[1])
        if i == sdsshort:
            L.append(j)
            print 'JUHU ',i,j
    sDict[i] = StatRot()
            #r = os.path.join('/',*t.split('/')[:-1])
            #if not fnmatch.fnmatch(os.path.basename(r), '*HZ.D'):
            #    print r

''
'''
'''
for root,dirs,files in os.walk(os.path.join(evpath,'data')):
        L = []
        for i in files:
            t = os.path.join(root,i)
            r = os.path.join('/',*t.split('/')[:-1])

            if not fnmatch.fnmatch(os.path.basename(r), '*HZ.D'):
                sdsname = i[:-11]
                L.append(t)
                for station in Meta:
                    if station.getcmpName() == sdsname:
                        print 'JUP',' SDSNAME ',sdsname,' station ',station
                        s_lat = float(station.lat)
                        s_lon = float(station.lon)
                        dist,azi,bazi = util.geodetics.gps2DistAzimuth(s_lat,s_lon,o_lat,o_lon)
                        print station,dist,azi,bazi,t
                        sDict[station] = StatRot(station.getName(),bazi,L)


for i in sDict.iterkeys():
    print i
'''

#rotate.rotate_NE_RT(n, e, ba)
