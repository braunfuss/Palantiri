
import os
import sys

#   add local directories to import path
                   
sys.path.append ('Common/')

#   import from Common

import Basic
import Globals
import NewVersion


def Usage () : return 'arraytool.py process event_name'

class Intern(object) :

    def __init__(self) : dummy = 0

    def error (self, text) :
        print '\nError : ' + text + '\n'
        print 'Usage : ' + Usage()
        sys.exit ('\n*** Program aborted ***')

    def checkProgramParameter (self, nMinParams, nMaxParams) :

        eventName = sys.argv[2]

        NewVersion.check ()                                        # Check if software version(s) ok

        if len (sys.argv) < nMinParams : self.error ('event name missing')
        if len (sys.argv) > nMaxParams : self.error ('Too many parameters')

        if not Globals.checkEventDirParameter (eventName):          # Exists event directory ?
           
           s =  'Invalid parameter - <' + eventName +'>'
           s += '\n        '
           s += 'Cannot find directory ' + Globals.EventDir()
           self.error (s)

#end class

# -------------------------------------------------------------------------------------------------

def start (config) :

    intern = Intern()

    if sys.argv[1] == 'process':
       
       #  python arraytool.py (0)  process (1) <event dir> (2)  
         
       intern.checkProgramParameter (3,3)

       path    = Globals.EventDir()
       at      = os.path.join (os.getcwd(),'Process', 'main.py')      # python script
       workDir = [path, 'tmp2', 'process']          # ???
       workDir = ['tmpProcess']
       cmd     = sys.executable + ' ' + at + ' -f ' + path

    else : 
       return False  

    Basic.changeDirectory (workDir)         # create and/or change working directory 
   #Basic.removeFiles     ('.')             # ... and empty it

    os.system (cmd)
    return True


