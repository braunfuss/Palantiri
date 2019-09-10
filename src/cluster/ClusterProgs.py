import os
import sys
from palantiri.common import Basic
from palantiri.common import  Globals

def Usage():
    return 'arraytool.py cluster event_name [automatic array clustering and print arrayconfig]'


class Intern(object):

    def __init__(self): dummy = 0

    def error(self, text):
        print('\nError: ' + text + '\n')
        print('Usage: ' + Usage())
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


def start(config):

    intern = Intern()

    if sys.argv[1] == 'cluster':
       intern.checkProgramParameter(3, 4)

       workDir = [Globals.EventDir(), 'tmp2', 'cluster']                   # ???
       workDir = ['cluster']
       cmd = "palantiri_cluster" + ' -f ' + Globals.EventDir()

    else:
       return False

    Basic.changeDirectory(workDir)         # create working directory

    os.system(cmd)
    return True
