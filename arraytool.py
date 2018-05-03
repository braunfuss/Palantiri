
import os
import fnmatch
import sys
from   optparse import OptionParser
import logging
import imp
import obspy.core

# add local directories to import path                      #hs+

sys.path.append ('Waveform/')
sys.path.append ('Process/')
sys.path.append ('Cluster/')
sys.path.append ('Common/')

import WaveformProgs
import ProcessProgs
import ClusterProgs
import CommonProgs
                                                            #hs-

config = imp.load_source ('Config',os.path.join (os.getcwd(),'tools','config.py'))

logger = logging.getLogger(sys.argv[0])
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

usage = ''' %prog [options] <command> [args]
%prog list              [list all events from event folder]
%prog search            [search for events in online catalog parameter from global.conf]
%prog create  eventID   [eventID from online catalog]
%prog getstations eventname
%prog getdata eventname
%prog getmeta eventname
%prog plotstations eventname
%prog cluster eventname     [automatic array clustering and print arrayconfig]
%prog process eventname
%prog plotresult eventname
'''

parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()

def folderContent(p):
    '''
    method to lookup necessary config files for event processing in the event folders

    return list of eventname if config and origin file are existing
    '''
    L = []
    for root,dirs,files in os.walk(p):
           flags=0

           for i in files:
                 if fnmatch.fnmatch (i,'*.config'):  flags += 1
                 if fnmatch.fnmatch (i,'*.origin'):  flags += 1

           if flags == 2:
               name = root.split('/')
               L.append(name[-1:])
    return L


def listEvents():
    '''
    method to list events in the event folder
    only list event if config and origin files are exisiting in event folder

    return list of events and print them
    '''

    for item in os.listdir(os.path.join(os.getcwd(),'events')):
        print item

# -------------------------------------------------------------------------------------------------

def parseArguments (args):
    '''
    parse arguments of main and entry script for arraytool

    do what you gave as commandline argument to main script
    '''

    dir = 'tools'

    #hs+
    if WaveformProgs.start (config) :  return       #  sys.argv[1] : 'getstations' 'getdata'  'getmeta'
    if ProcessProgs.start  (config) :  return       #                'process'
    if ClusterProgs.start  (config) :  return       #                'cluster'
    if CommonProgs.start   () :        return       #                'new_version  (for internal use)
    #hs-

    if sys.argv[1] == 'list':
        listEvents()

    elif sys.argv[1] == 'process':            # nicht mehr benutzt

        path = os.path.join (os.getcwd(),'events',sys.argv[2])
        at   = os.path.join (os.getcwd(), dir,      'main.py')
        cmd  = sys.executable+' '+at+' -f '+ path

        os.chdir  (os.path.join(os.getcwd(),"tools"))
        os.system (cmd)

    elif sys.argv[1] == 'energy':

        path = os.path.join (os.getcwd(),'events',sys.argv[2])
        at   = os.path.join (os.getcwd(), dir,    'energy.py')
        cmd  = sys.executable+' '+at+' -f '+ path

        os.chdir  (os.path.join (os.getcwd(),"tools"))
        os.system (cmd)

    elif sys.argv[1] == 'create':

        at = os.path.join (os.getcwd(),dir,'create.py')
        t  = ''

        for i in sys.argv[2:] : t +=i+'_'

        cmd = sys.executable+' '+at+' '+t
        os.chdir (os.path.join (os.getcwd(),"tools"))
        os.system(cmd)
        '''

    elif sys.argv[1] == 'rotate':

        at = os.path.join(os.getcwd(),dir,'rotatecomp.py')
        path = os.path.join(os.getcwd(),'events',sys.argv[2])
        print path
        cmd = sys.executable+' '+at+' -p '+ path

        os.chdir (os.path.join(os.getcwd(),"tools"))
        os.system(cmd)
        '''

    elif sys.argv[1] == 'cluster':

        at   = os.path.join(os.getcwd(), dir,     'callcluster.py')
        path = os.path.join(os.getcwd(),'events',sys.argv[2])
        cmd  = sys.executable+' '+at+' -f '+ path

        os.chdir (os.path.join (os.getcwd(),"tools"))
        os.system(cmd)

    elif sys.argv[1] == 'plotstations':

        at   = os.path.join(os.getcwd(),dir,'makeplot.py')
        path = os.path.join(os.getcwd(),'events',sys.argv[2])
        cmd  = sys.executable+' '+at+' -f '+ path

        os.chdir (os.path.join(os.getcwd(),"tools"))
        os.system(cmd)

    elif sys.argv[1] == 'pyrocko_download':

        at   = os.path.join(os.getcwd(),'Waveform','pyrocko_down.py')
        path = os.path.join(os.getcwd(),'events',sys.argv[2])
        cmd  = sys.executable+' '+at+' -f '+ path

        os.chdir (os.path.join(os.getcwd(),"Waveform"))
        os.system(cmd)


    elif sys.argv[1] == 'search':

        at  = os.path.join(os.getcwd(),dir,'eventsearch.py')
        cmd = sys.executable+' '+ at

        os.system (cmd)

    elif sys.argv[1] == 'plotresult':

        at   = os.path.join(os.getcwd(),dir,'plotresult.py')
        path = os.path.join(os.getcwd(),'events',sys.argv[2])
        cmd  = sys.executable+' '+at+' -f '+ path

        os.chdir(os.path.join(os.getcwd(),"tools"))
        os.system(cmd)


    else:
        logger.info('\033[31m Option not available \033[0m')

# -------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    if len(sys.argv) > 1:
        parseArguments(sys.argv)
    else:
        cmd = 'python '+sys.argv[0]+' --help'
        os.system(cmd)
