#
#      Implements own global data
#
import os
import sys
import getpass
import platform

WINDOWS =(platform.system() == 'Windows')

#      Imports from Common

import Basic
import Logfile
import ConfigFile

#      Constants

#      Global variables

isClient   = False                        # is client process ?
isDebug= False                        # debug mode ?

ConfigDict  = None                        # Configuration file as dictionary   
ProtFileDir = None                        # directory for protocol files

EVENTS = 'events'

_eventDir = None

# -------------------------------------------------------------------------------------------------

def setEventDir(s):
    global _eventDir

    _eventDir = s
    return s


def EventDir():                          # event directory
    global _eventDir

    if _eventDir == None:
       #print 'Globals: event dir not set'
       #assert False

       s = os.path.join(os.getcwd(), EVENTS)
       n = len(sys.argv)

       if   n < 3 : _eventDir = s
#      else:        _eventDir = os.path.join(s, sys.argv[n-1])
       else:        _eventDir = os.path.join(s, sys.argv[2])
    #endif

    return _eventDir

def TempFileName(name):

    assert ProtFileDir != None
    Basic.createDirectory(ProtFileDir)
    return os.path.join(ProtFileDir,name)

def KeyfileFolder():
    return os.path.join(EventDir(), ConfigDict ['keyfilefolder'])

# -------------------------------------------------------------------------------------------------

def _error(text):
    print '\nError: ' + text + '\n'
    sys.exit(1)


def checkNrOfParameter(nMin, nMax):

    if len(sys.argv) < nMin: _error('event name missing')
    if len(sys.argv) > nMax: _error('Too many parameters')


# Exists event directory ?

def checkEventDirParameter(param):

    s1  = os.path.basename(param)
    dir = os.path.join(os.getcwd(), EVENTS, s1)
    #print 'dir=', dir

    return os.path.isdir(dir)
 
# -------------------------------------------------------------------------------------------------

def init(configFileName = None):

    global EventDir, ProtFileDir, ConfigDict, isDebug

    ProtFileDir = os.path.join(EventDir(), 'tmp1')

    if True:   # not isClient:
       ConfigDict  = ConfigFile.readGlobalConf(configFileName)

       if ConfigDict == None:  return False

    #obj = ConfigFile.GlobCfg()
    #isDebug = obj.Bool('debug', '0')
    key = 'DEBUG_FLAG'

    if not os.environ.has_key(key): isDebug = False
    else:                            isDebug =(os.environ [key].strip() == '1')

    if not isClient:
       if isDebug: Logfile.add('Debugging is on')

    return True

