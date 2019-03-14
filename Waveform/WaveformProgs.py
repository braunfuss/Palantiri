import os
import sys
import obspy.core

# add local directories to import path

sys.path.append('Common/')

#      Import from Common

import Globals


def Usage(action): return 'arraytool.py ' + action + ' event_name'

class Intern(object):

    def __init__(self, action):  self.action = action

    def error(self, text):

        print('\nError: ' + text + '\n')
        print('Usage: ' + Usage(self.action))
        sys.exit('\n*** Program aborted ***')

    def checkProgramParameter(self, nMinParams, nMaxParams):

        eventName = sys.argv[2]


        if len(sys.argv) < nMinParams: self.error('event name missing')
        if len(sys.argv) > nMaxParams: self.error('Too many parameters')

        if not Globals.checkEventDirParameter(eventName):          # Exists event directory ?

           s =  'Invalid parameter - <' + eventName +'>'
           s += '\n        '
           s += 'Cannot find directory ' + Globals.EventDir()
           self.error(s)

#end class

# -------------------------------------------------------------------------------------------------

def buildCommand(config):

    dir= 'Waveform'                                  # directory of local python scripts
    action = sys.argv[1]
    intern = Intern(action)

    if action == 'getstations':

        #  python arraytool.py(0)  getmeta(1) <event dir>(2)

        intern.checkProgramParameter(3,3)

        at   = os.path.join(os.getcwd(), dir, 'getStationList.py')
        cmd  = sys.executable + ' ' + at + ' -f ' + Globals.EventDir()

    elif action == 'getdata':

        #  python arraytool.py(0)  getdata(1) <event dir>(2)  [network(3)]

        intern.checkProgramParameter(3,4)

        at   = os.path.join(os.getcwd(), dir, 'getStationWaveformData.py')
        cmd  = sys.executable + ' ' + at + ' -f ' + Globals.EventDir()

        if len(sys.argv) > 3:  cmd +=(' -n ' + sys.argv[3])                # specific network

    elif action == 'getmeta':

        #  python arraytool.py(0)  getmeta(1) <event dir>(2)  [network(3)]

        intern.checkProgramParameter(3,4)

        path   = Globals.EventDir()                                          # event directory
        at = os.path.join(os.getcwd(), dir, 'ev_meta_mt4.py')
        C  = config.Config(path)
        Origin = C.parseConfig('origin')

        d   = obspy.core.utcdatetime.UTCDateTime(Origin['time'])
        jd  = "%03d" % d.julday
        cmd =('%s %s -p %s -y %s -d %s' )%(sys.executable,at, path, str(d.year), str(jd))

        if len(sys.argv) > 3:  cmd +=(' -n ' + sys.argv[3])               # specific network

    else:
       return None

    return cmd

# -------------------------------------------------------------------------------------------------

def start(config):

    cmd = buildCommand(config)

    if cmd == None: return False

    os.chdir(os.path.join(os.getcwd(), "tools"))                 # ??? Directory aendern
    os.system(cmd)

    return True
