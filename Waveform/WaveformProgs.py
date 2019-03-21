import os
import sys
sys.path.append('Common/')
import Globals


def Usage(action):
    return 'arraytool.py ' + action + ' event_name'


class Intern(object):

    def __init__(self, action):  self.action = action

    def error(self, text):

        print('\nError: ' + text + '\n')
        print('Usage: ' + Usage(self.action))
        sys.exit('\n*** Program aborted ***')

    def checkProgramParameter(self, nMinParams, nMaxParams):

        eventName = sys.argv[2]

        if len(sys.argv) < nMinParams:
            self.error('event name missing')
        if len(sys.argv) > nMaxParams:
            self.error('Too many parameters')

        if not Globals.checkEventDirParameter(eventName):

            s = 'Invalid parameter - <' + eventName +'>'
            s += '\n        '
            s += 'Cannot find directory ' + Globals.EventDir()
            self.error(s)


def buildCommand(config):

    return None


def start(config):

    cmd = buildCommand(config)

    if cmd is None:
        return False

    os.chdir(os.path.join(os.getcwd(), "tools"))
    os.system(cmd)

    return True
