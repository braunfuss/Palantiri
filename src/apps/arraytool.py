import os
import fnmatch
import sys
from optparse import OptionParser
import logging
import imp
import obspy.core
from palantiri.process import ProcessProgs
from palantiri.cluster import ClusterProgs
from palantiri.common import CommonProgs
from palantiri.tools import config as config

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
%prog pyrocko_down eventname
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
    for root,dirs ,files in os.walk(p):
        flags = 0

        for i in files:
            if fnmatch.fnmatch(i, '*.config'):
                flags += 1
            if fnmatch.fnmatch(i, '*.origin'):
                flags += 1

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
        print(item)



def parseArguments(args):
    '''
    parse arguments of main and entry script for arraytool

    do what you gave as commandline argument to main script
    '''

    dir = 'tools'

    if ProcessProgs.start(config):
        return
    if ClusterProgs.start(config):
        return


    if sys.argv[1] == 'list':
        listEvents()


    elif sys.argv[1] == 'create':

        at = "palantiri_create"
        t  = ''

        for i in sys.argv[2:] : t +=i+'_'

        cmd = at+' '+t
        os.system(cmd)

    elif sys.argv[1] == 'cluster':

        at = "palantiri_cluster"
        cmd  = at+' -f '+ path
        os.system(cmd)

    elif sys.argv[1] == 'pyrocko_download':

        at   = os.path.join(os.getcwd(),'Waveform','pyrocko_down.py')
        path = os.path.join(os.getcwd(),'events',sys.argv[2])
        cmd  = sys.executable+' '+at+' -f '+ path

        os.chdir(os.path.join(os.getcwd(),"Waveform"))
        os.system(cmd)

    elif sys.argv[1] == 'search':
        os.system("palantiri_eventsearch")

    else:
        logger.info('\033[31m Option not available \033[0m')

def main():
    if len(sys.argv) > 1:
        parseArguments(sys.argv)
    else:
        cmd = 'python '+sys.argv[0]+' --help'
        os.system(cmd)
