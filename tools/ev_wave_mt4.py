from optparse import OptionParser
from obspy.core import read,Stream
import os
import fnmatch
import fileinput
import re
import obspy.core
import logging
import multiprocessing
import time
import subprocess
import sys
import config
from ConfigParser import SafeConfigParser
from obspy.arclink.client import Client
import urllib
import urllib2
from obspy.core.util import NamedTemporaryFile

logger = logging.getLogger(sys.argv[0])
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

parser = OptionParser(usage="%prog -t 2009-12-31T12:23:34 -d 5 -m SDS -k key -s all/acq")
parser.add_option("-t","--time", type="string", dest="time", help="time")
parser.add_option("-d","--duration", type="string", dest="duration", help="duration in min")
parser.add_option("-m","--sdsfolder", type="string", dest="sdsfolder", help="sdsfolder")
parser.add_option("-k","--keyfolder", type="string", dest="keyfolder", help="keyfolder")
parser.add_option("-s","--station", type="string", dest="stationversion", help="stationversion")
parser.add_option("-f","--evpath", type="string", dest="eventpath", help="eventpath")

(options, args) = parser.parse_args()

#_key_search_folder=options.keyfolder
#_mseed_search_folder=options.sdsfolder

def globalConf():
    cDict = {}
    parser = SafeConfigParser()
    parser.read('../global.conf')
    for section_name in parser.sections():
        for name, value in parser.items(section_name):
            cDict[name]=value
    return cDict

def make_time(begin_time,duration):

    b_time_1 = obspy.core.utcdatetime.UTCDateTime(begin_time,iso8601=True) 
    geofon_begin = b_time_1.formatArcLink()

    e_time = b_time_1 + float(duration) * 60
    geofon_end = e_time.formatArcLink()
    
    b_time_iris = obspy.core.utcdatetime.UTCDateTime(begin_time,iso8601=True) 
    iris_begin = b_time_iris.formatIRISWebService()

    e_time_iris = b_time_iris + float(duration) * 60
    iris_end = e_time_iris.formatIRISWebService()
    
    dict_time = {'g_begin':geofon_begin,'g_end':geofon_end,
                 'i_begin':iris_begin,'i_end':iris_end,
                 'obspy_begin':b_time_1,'obspy_end':e_time}

    return dict_time

def keyfolder2List(p):

    logger.info('\033[31m Parsing KEYFILE Structure \033[0m \n')
    
    stationfilespath=p
    L = []

    pattern = 'station_' + "[A-Za-z]"
    content = os.listdir(stationfilespath)

    for i in content:
        if re.match(pattern, i):
            name = i[8:]
            L.append(name)

    L_sauber = list(set(L))
    L_sauber.sort()

    K = []
    for i in L_sauber:
        if options.stationversion == 'acq':
            f = 'station_'+i
            path = os.path.join(stationfilespath,f)
            for line in fileinput.input(path):
                if "PACKAGES" in line:
                    if fnmatch.fnmatch(line, '*acquisition*'):
                        K.append(i)
        if options.stationversion == 'all':
            K.append(i)
        
        if options.stationversion == '':
            print 'set station method'
        
    logger.info('\033[31m FOUND: --> %d KEYFILES \033[0m\n' % (len(K)))
    
    return K

def sdsList(p):

    L = []

    logger.info('\033[31m Parsing SDS Structure \033[0m\n')
    for root,dirs,files in os.walk(p):
        for i in files:
            line = str.split(i,'.')
            val= line[0]+'_'+line[1]
            L.append(val)

    L_sauber = list(set(L))

    logger.info('\033[31m FOUND ---> %d STATIONS \033[0m\n' % (len(L_sauber)))
    
    return L_sauber

def cmpSdsKey(K,S):

    logger.info('\033[31m Parsing for Missing Station \033[0m\n')
    t = list(set(K).difference(set(S)))
    t = set(t)
    logger.info('\033[31m %d STATIONS are missing \033[0m' % (len(t)))

    return t


def to_sds(output):
    
    cmd = 'scart -I '+output+' '+_mseed_search_folder
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()
    logger.info('\033[31munpacking %s to SDS\n\n \033[0m' % (cmd))

def create_dir(dir):
  if os.access(dir, os.W_OK):
    return True
  try:
    os.makedirs(dir)
    return True
  except:
    return False


def multiplex(filename):

    st = read(filename)
    for trace in st:
        path = os.path.join(_mseed_search_folder,str(trace.stats.starttime.year),trace.stats.network,trace.stats.station,('%s.D')%(trace.stats.channel))
        create_dir(path)

        diff = trace.stats.endtime.day - trace.stats.starttime.day
        rangediff = diff+1
        if diff == 0:
            d = obspy.core.utcdatetime.UTCDateTime(trace.stats.starttime)
            jd = str("%03d" % d.julday)
            finalname = ('%s.%s.%s.%s.D.%s.%s')%(trace.stats.network,trace.stats.station,trace.stats.location,trace.stats.channel,trace.stats.starttime.year,jd)
            filepath = os.path.join(path,finalname)
            trace.write(filepath,format='MSEED',reclen=512)
        elif diff >= 1:
            for i in xrange(rangediff):
                mult = 60*60*24*i
                if i== 0:
                    s1 = trace.stats.starttime
                    e1 = obspy.core.utcdatetime.UTCDateTime(year=s1.year, month=s1.month, day=s1.day,hour = 23, minute=59,second=59, microsecond=999999)
                else:
                    s1 = trace.stats.starttime+mult
                    s1 = obspy.core.utcdatetime.UTCDateTime(year=s1.year, month=s1.month, day=s1.day,hour = 0, minute=0, microsecond=00)
                    e1 = obspy.core.utcdatetime.UTCDateTime(year=s1.year, month=s1.month, day=s1.day,hour = 23, minute=59,second=59, microsecond=999999)
                if i == diff:
                    e1 = trace.stats.endtime
                d = obspy.core.utcdatetime.UTCDateTime(s1)
                jd = str("%03d" % d.julday)
                finalname = ('%s.%s.%s.%s.D.%s.%s')%(trace.stats.network,trace.stats.station,trace.stats.location,trace.stats.channel,trace.stats.starttime.year,jd)
                filepath = os.path.join(path,finalname)
                tr = trace.slice(s1, e1)
                tr.write(filepath,format='MSEED',reclen=512)
                print s1,' to ',e1,jd,finalname,filepath


def proof_file_v2(filename,component):
    
    size=0
    streamList = []
    
    if (os.path.isfile(filename) and os.path.getsize(filename) != 0):
          
          try:
              st = read(filename)
              for i in st:
                  streamList.append(i.stats.channel)
              if len(set(streamList)) == len(component):
                  #to_sds(filename)
                  multiplex(filename)
                  size = os.path.getsize(filename)
              else:
                  size = 0
          except:
               size = 0
    
    return size

def make_arcObspyrequest(station,begin,end,channel,loc):
    
    logger.info('\033[31m Obspy ARCLINKREQUEST \033[0m\n' % ())
    size = 0
    pid = os.getpid()

    rstation = station.split('_')
    
    if loc == '--':
        loc = ''
    client = Client(user=_mail,dcid_keys={'GFZ': '5]p[x#mJ'},timeout=10)

    st = obspy.core.stream.Stream()
    for i in channel:
        st += client.getWaveform(rstation[0], rstation[1], loc, i, begin, end)
    
    
    output = str(pid)+'-'+station+'-ARCLINK.mseed'
    st.write(output,format='MSEED',reclen=512)
    size = proof_file_v2(output,channel)
    
    try:
        os.remove(output)
    except:
        print ''
    
    return size


def make_arcrequest(station,begin,end,channel,loc):

    logger.info('\033[31m ARCLINKREQUEST \033[0m\n' % ())
    rstation = station.replace('_',' ')
    pid = os.getpid()
    reqfile = 'ARC-'+str(pid)+'-'+station+'-'+loc+'.txt'
    fname = os.path.join(os.getcwd(),reqfile)
    fobj = open(fname,'w')
    if loc == '--':
        loc = ''
    for a in channel:
            content = str(begin)+' '+str(end)+' '+rstation+' '+a+' '+loc+'\n'
            logger.info('%s %s %s %s %s' % (rstation, loc, a, str(begin) ,str(end)))
            fobj.write(content)
    fobj.close()
    
    output = str(pid)+'-'+station+'-ARCLINK.mseed'
    cmd=''
    if fnmatch.fnmatch(station, 'MN_*'):
        cmd = 'arclink_fetch -qqq -p -a webdc.eu:18002 -w dcidpasswords.txt -k mseed -u '+_mail+' --timeout=10 -o '+output+' '+fname
    else:
        cmd = 'arclink_fetch -qqq -a webdc.eu:18002 -w dcidpasswords.txt -k mseed -u '+_mail+' --timeout=10 -o '+output+' '+fname

    #p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #p.wait()
    
    os.system(cmd)
    size = proof_file_v2(output,channel)
    
    try:
        os.remove(fname)
        os.remove(output)
    except:
        print ''
    
    return size


def make_irisrequest(station,begin,end,channel,loc):
    
    logger.info('\033[31m IRISREQUEST \033[0m\n' % ())
    
    net,sta = station.split('_')
    url='http://service.iris.edu/fdsnws/dataselect/1/query?'
    st = Stream()
    for i in channel:
        parameter = urllib.urlencode({
                 'net': net,
                 'sta': sta,
                 'loc': loc,
                 'cha': i,
                 'starttime': begin,
                 'endtime': end,
                 'nodata': '404',
        })
        u = ('%s%s')%(url,parameter)
        data = urllib.urlopen(u).read()
        tf = NamedTemporaryFile()
        tf.write(data)
        tf.seek(0)
        st += read(tf.name, 'MSEED')
        tf.close()
        
    print st
    output = ('%s.%s-IRIS.mseed')%(net,sta)
    st.write(output,format='MSEED')

    size = proof_file_v2(output,channel)
    #print 'SIZE: ----> ',size
    try:
        #os.remove(fname)
        os.remove(output)
    except:
        print ''

    return size
    

def write_statistic(year,day):

    K = keyfolder2List(_key_search_folder)
    S = sdsList(_mseed_search_folder)
    MISS = cmpSdsKey(K,S)

    fobjm = open(os.path.join(options.eventpath,'missing-station-'+str(year)+'-'+str(day)+'.dat'),'w')
    fobjm.write('\n'.join(MISS))
    fobjm.close()
    
    fobja = open(os.path.join(options.eventpath,'available-station-'+str(year)+'-'+str(day)+'.dat'),'w')
    fobja.write('\n'.join(S))
    fobja.close()

def make_meta(stime):
    
    d = obspy.core.utcdatetime.UTCDateTime(stime)
    jd = "%03d" % d.julday
    cmd='python ev_meta_mt4.py -p '+_mseed_search_folder+' -y '+str(d.year)+' -d '+str(jd)
    write_statistic(d.year,jd)
    print cmd
    
    #p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #p.wait()
    #os.system(cmd)

def stationRequest(station,timeDict,counter):
    
    #signal.signal(signal.SIGINT, signal_handler)
    for i in _loc:
        for j in _channel:
            try:
                logger.info('\033[31m COUNTER: %d %s %s %s %s %s \033[0m\n' % (counter,station,i,j,timeDict['i_begin'],timeDict['i_end']))
                size = make_irisrequest(station,timeDict['i_begin'],timeDict['i_end'],j,i)
                if size != 0:
                    logger.info('\033[31m DATA FROM IRIS FOR STATION % LOC: %s CHANNEL: %s \033[0m\n' % (station,i,j))
                    return
                else:
                    #size = make_arcrequest(station,timeDict['g_begin'],timeDict['g_end'],j,i)
                    #size = make_arcObspyrequest(station,timeDict['obspy_begin'],timeDict['obspy_end'],j,i)
                    pass
                    if size != 0:
                        logger.info('\033[31m DATA FROM WEBDC FOR STATION % LOC: %s CHANNEL: %s \033[0m\n' % (station,i,j))
                        return
                    else:
                        logger.info('\033[31m NO DATA FROM WEBDC OR IRIS FOR STATION %s AND FINISH  \033[0m\n' % (station))
            except:
                continue



def run():
    
    time_d = make_time(options.time,options.duration)
    K = keyfolder2List(options.keyfolder)
    
    T = []
    S = sdsList(options.sdsfolder)
    MISS = cmpSdsKey(K,S)

    for i in MISS:
        T.append(i)
    
    po = multiprocessing.Pool(1)
    for i in xrange(len(MISS)):
        po.apply_async(stationRequest,(T[i],time_d,i))
    po.close()
    po.join()
    
    
if __name__ == "__main__":

    channel1 = ['BHE','BHN','BHZ']
    channel2 = ['BH1','BH2','BHZ']
    channel3 = ['HHE','HHN','HHZ']
    _channel = [channel1,channel2,channel3]
    _loc = ['--','10','00','11']
    
    C = config.Config(options.eventpath)
    Origin = C.parseConfig('origin')
    Conf = globalConf()
    
    _mail = Conf['mail']
    #print options.eventpath,Conf
    
    options.time=Origin['time']
    options.duration=int(Conf['duration'])
    options.sdsfolder=os.path.join(options.eventpath,'data')
    options.keyfolder=os.path.join(options.eventpath,Conf['keyfilefolder'])
    options.stationversion = 'all'
    
    _key_search_folder =options.keyfolder
    _mseed_search_folder = options.sdsfolder 

    t1 = time.time()
    run()
    t2 = time.time()
    
    print '%s took %0.3f s' % ('CALC :', (t2-t1))
    
    make_meta(options.time)
