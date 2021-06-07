import os
from configparser import SafeConfigParser
import sys
import shutil
import logging
from urllib.request import urlopen
import dateutil.parser
from obspy.core.utcdatetime import UTCDateTime

logger = logging.getLogger('ARRAY-MP')
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s -\
                               %(levelname)s - %(message)s")
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)


def init():
    '''
    method to parse global conf

    return Configdictionary
    '''
    confname = os.path.join(os.getcwd(), '..', 'global.conf')

    cDict = {}
    parser = SafeConfigParser()
    parser.read(confname)
    for section_name in parser.sections():
        for name, value in parser.items(section_name):
            cDict[name] = value

    logger.info('\033[31m Global Configuration %s \033[0m \n' % (cDict))
    return cDict


def parseEvent(eventID):

    eventID = eventID[1].replace('_', '')

    url = 'http://service.iris.edu/fdsnws/event/1/query?eventid=%s&format=text' % (eventID)
    data = urlopen(url).read()
    data = data.decode('utf_8')
    data = data.split('\n')
    dl = data[1:]
    i = dl[0]
    i = i.split('|')
    time = i[1].replace(':', '-').strip()
    name = i[12].replace(' ', '-').strip()
    eventname = ('%s_%s') % (name, time)

    return eventname


def createWorkingDirectory(args):
    '''
    method to create working directory of event to be processed

    return folder path and event_id
    '''
    foldername = parseEvent(args)
    absfolder = os.path.join(os.getcwd(), 'events', foldername)

    if os.access(absfolder, os.F_OK) is False:
        os.makedirs(absfolder)
        logger.info('\033[31m WORKING FOLDER CREATED \033[0m \n')

    logger.info('\033[31m Folder: %s  EventID: %s \033[0m \n' % (absfolder,
                foldername))

    return absfolder, foldername


def writeOriginFile(path, ev_id):
    '''
    method to write origin(event) file in the event directory to be processed

    return origin time of event
    '''
    fname = os.path.basename(path)+'.origin'
    fobj = open(os.path.join(path, fname), 'w')
    fobj.write('[origin]\n\n')

    eventID = ev_id[1].replace('_', '')

    url = 'http://service.iris.edu/fdsnws/event/1/query?eventid=%s&format=text' % (eventID)
    data = urlopen(url).read()
    data = data.decode('utf_8')
    data = data.split('\n')
    dl = data[1:]
    i = dl[0]
    i = i.split('|')
    time = str(dateutil.parser.parse(i[1]))[:19]
    fobj.write('region = %s\n' % i[12].strip())
    fobj.write('lat = %s\n' % i[2])
    fobj.write('lon = %s\n' % i[3])
    fobj.write('depth = %s\n' % i[4])
    fobj.write('time = %s\n' % time)
    fobj.write('strike = -999\n')
    fobj.write('dip = -999\n')
    fobj.write('rake = -999\n')

    fobj.close()

    return time


def writeSynFile(path, ev_id):
    '''
    method to write synthetic input(event) file in the event directory to be
    processed

    '''
    fname = os.path.basename(path)+'.syn'
    fobj = open(os.path.join(path, fname), 'w')
    fobj.write('[synthetic parameter]\n\n')

    eventID = ev_id[1].replace('_', '')

    url = 'http://service.iris.edu/fdsnws/event/1/query?eventid=%s&format=text' % (eventID)
    data = urlopen(url).read()
    data = data.decode('utf_8')
    data = data.split('\n')
    dl = data[1:]
    i = dl[0]
    i = i.split('|')
    time = i[1]
    fobj.write('region = %s\n' % i[12].strip())
    fobj.write('nsources = 1\n')
    fobj.write('lat_0 = %s\n' % i[2])
    fobj.write('lon_0 = %s\n' % i[3])
    fobj.write('depth_0 = %s\n' % i[4])
    fobj.write('time_0 = %sZ\n' % i[1])
    fobj.write('strike_0 = -999\n')
    fobj.write('dip_0 = -999\n')
    fobj.write('rake_0 = -999\n')
    fobj.write('width_0 = -999\n')
    fobj.write('length_0 = -999\n')
    fobj.write('slip_0 = -999\n')
    fobj.write('nucleation_x_0 = 0\n')
    fobj.write('nucleation_y_0 = 0\n')
    fobj.write('store = store_id\n')
    fobj.write('store_superdirs = dir of store\n')
    fobj.write('use_specific_stf = 0\n')
    fobj.write('stf = stf=gf.HalfSinusoidSTF()\n')
    fobj.write('source = RectangularSource\n')

    fobj.close()

    return time


def copyConfigSkeleton(evfolder):
    '''
    method to copy the example config from skeleton directory to the
    event directory
    '''
    logger.info('\033[31m Copy example.config to %s \033[0m \n' % (evfolder))

    dstfile = os.path.split(evfolder)[1]+'.config'
    import palantiri
    path = palantiri.__path__
    src = os.path.join(path[0], 'skeleton', 'example.config')
    dst = os.path.join(evfolder, dstfile)
    logger.info('\033[31m Created event directory \
                %s \033[0m \n' % (evfolder.split('/')[-1]))
    shutil.copy(src, dst)

    event = evfolder.split('/')[-1]
    eventname = event.split('_')[0]
    time = event.split('_')[-1]
    time = UTCDateTime(time)
    time = time.format_iris_web_service()
    time_p1 = time[0:10]
    time_p2 = time[11:-1]
    time_pf = time_p1 + " " + time_p2

    logger.info('\033[31mNEXT PROCESSING STEP (e.g.): \n\n   \
                 palantiri_down {evdirectory} "%s" 10352.323104588522 0.1 10 --radius-min=1110 %s \n\n\033[0m'.format(evdirectory=str(event.strip('[]'))) %(time_pf, eventname))


def main():

    if len(sys.argv) == 2:
        absf, evid = createWorkingDirectory(sys.argv)
        writeSynFile(absf, sys.argv)
        writeOriginFile(absf, sys.argv)
        copyConfigSkeleton(absf)
    else:
        logger.info('\033[31m Nothing to do %s \033[0m \n')
