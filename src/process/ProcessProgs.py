import os
import sys
from palantiri.common import Basic
from palantiri.common import Globals


def Usage():
    return 'arraytool.py process event_name'


class Intern(object):

    def __init__(self):
        dummy = 0

    def error(self, text):
        print('\nError: ' + text + '\n')
        print('Usage: ' + Usage())
        sys.exit('\n*** Program aborted ***')

    def checkProgramParameter(self, nMinParams, nMaxParams):

        eventName = sys.argv[2]

        if len(sys.argv) < nMinParams:
            self.error('event name missing')
        if len(sys.argv) > nMaxParams:
            self.error('Too many parameters')

        if not Globals.checkEventDirParameter(eventName):

            s = 'Invalid parameter - <' + eventName + '>'
            s += '\n        '
            s += 'Cannot find directory ' + Globals.EventDir()
            self.error(s)


def start(config):

    intern = Intern()

    if sys.argv[1] == 'process':
        intern.checkProgramParameter(3, 4)

        path = Globals.EventDir()
        path_emp = Globals.EventDir_emp()
        try:
            path_emp = Globals.EventDir_emp()
            at = os.path.join(os.getcwd(), 'Process', 'main.py')
            workDir = [path, 'tmp2', 'process']
            workDir = ['tmpProcess']
            cmd = "palantiri_process" + '  -f ' + path + ' -e ' + path_emp
        except Exception:
            at = os.path.join(os.getcwd(), 'Process', 'main.py')
            workDir = [path, 'tmp2', 'process']
            workDir = ['tmpProcess']
            cmd = "palantiri_process" + ' -f ' + path
    else:
        return False

    Basic.changeDirectory(workDir)
    os.system(cmd)
    return True
