import os
import sys
import getpass
import platform
import Basic
import Logfile
import ConfigFile

isClient = False
isDebug= False

ConfigDict  = None                        # Configuration file as dictionary
ProtFileDir = None                        # directory for protocol files

EVENTS = 'events'

_eventDir = None
_eventDir_emp = None


def setEventDir(s):
    global _eventDir

    _eventDir = s
    return s

def setEventDir_emp(s):
    global _eventDir_emp

    _eventDir_emp = s
    return s

def EventDir():                          # event directory
    global _eventDir

    s = os.path.join(os.getcwd(), EVENTS)
    n = len(sys.argv)

    eventDir = os.path.join(s, sys.argv[2])

    return eventDir


def EventDir_emp():                          # event directory
    global _eventDir_emp

    s = os.path.join(os.getcwd(), EVENTS)
    n = len(sys.argv)
    try:
        _eventDir_emp  = os.path.join(s, sys.argv[3])
    except IndexError:
        _eventDir_emp  = None
    return _eventDir_emp

def TempFileName(name):

    assert ProtFileDir != None
    Basic.createDirectory(ProtFileDir)
    return os.path.join(ProtFileDir,name)

def KeyfileFolder():
    return os.path.join(EventDir(), ConfigDict ['keyfilefolder'])

# -------------------------------------------------------------------------------------------------

def _error(text):
    print('\nError: ' + text + '\n')
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
    key = 'DEBUG_FLAG'

    if not isClient:
       if isDebug: Logfile.add('Debugging is on')

    return True
