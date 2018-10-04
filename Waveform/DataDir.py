
import os
import sys
import platform

WINDOWS =(platform.system() == 'Windows')

# add local directories to import path  
                  
sys.path.append('../Common/')

#import fnmatch

import  obspy.core.trace
 
#      Import from Common

import Globals                     # Own global data
import Basic                       # Own module with basic functions
import Logfile                     # Implements logfile
from   DataTypes import Station

DATA_DIR         = 'data'                   # root name of data directory(relativ to event dir)
FILE_NAME_FORMAT = '%s.%s.%s.%s.D.%s.%s'

# -------------------------------------------------------------------------------------------------

def filename(trace, day) :

    postfix = str("%03d" % day.julday)

    if type(trace) is obspy.core.trace.Trace :
       t        = trace.stats
       filename =(FILE_NAME_FORMAT) %(t.network, t.station, t.location, t.channel, 
                                        t.starttime.year, postfix)
    else :
       Logfile.exception('DataDir.filename', str(type(trace)))
       Logfile.abort('')

    #Logfile.add(filename)
    return filename

# -------------------------------------------------------------------------------------------------

def getFileNames(eventDir=None) :

    if eventDir == None :  eventDir = Globals.EventDir()

    names = []
    path  = os.path.join(eventDir, DATA_DIR)

    for root,dirs,files in os.walk(path):
        for s in files : names.append(s)
   
    #Logfile.addLines(names)
    return sorted(names)

# -------------------------------------------------------------------------------------------------

def getNetworks(eventDir=None) :

    files    = getFileNames(eventDir)
    networks = []

    for s in files :
        net = str.split(s, '.')[0]
        networks.append(net)

    networks = sorted(set(networks))
    #Logfile.addLines(networks)

    return networks

def isNetwork(network, eventDir=None) :

    assert network != None
    return network in getNetworks(eventDir)

# -------------------------------------------------------------------------------------------------



