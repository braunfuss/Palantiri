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
        sys.exit('\n*** Gandalf made Pippin drop the Palantiri by ***')

    def checkProgramParameter(self, nMinParams, nMaxParams):


        eventName = sys.argv[2]
        if len(sys.argv) < nMinParams:
            self.error('event name missing')
        if not Globals.checkEventDirParameter(eventName):
            s = 'Invalid parameter - <' + eventName + '>'
            s += '\n        '
            s += 'Cannot find directory ' + Globals.EventDir()
            self.error(s)


def start(config):
    intern = Intern()
    if sys.argv[1] == 'process':
        force = False
        for argv in sys.argv[1:]:
            if argv == '--force':
                force = True

        intern.checkProgramParameter(3, 4)

        path = Globals.EventDir()
        path_emp = Globals.EventDir_emp()
        try:
            path_emp = Globals.EventDir_emp()
            at = os.path.join(os.getcwd(), 'Process', 'main.py')
            workDir = [path, 'tmp2', 'process']
            workDir = ['tmpProcess']
            if force is False:
                cmd = "palantiri_process" + '  -f ' + path + ' -e ' + path_emp
            else:
                cmd = "palantiri_process" + '  -f ' + path + ' -e ' + path_emp + '--force'

        except Exception:
            at = os.path.join(os.getcwd(), 'Process', 'main.py')
            workDir = [path, 'tmp2', 'process']
            workDir = ['tmpProcess']
            if force is False:
                cmd = "palantiri_process" + ' -f ' + path
            else:
                cmd = "palantiri_process" + ' -f ' + path + '--force'

    else:
        return False

    Basic.changeDirectory(workDir)
    os.system(cmd)
    return True
