import os, sys, platform

WINDOWS = (platform.system() == 'Windows')

import Basic, Globals, Logfile                     # Import from Common

class MainObj (object) :

   def __init__ (self, externClass, version, runTimeLog=None, errorLog=None) :

       self.version= version
       self.runTimeLog = runTimeLog
       self.errorLog   = errorLog
       self.extern = externClass

   def run (self) :

       Logfile.init (self.runTimeLog, self.errorLog)
       Logfile.setStartMsg (self.version)

       if not self.extern.init() : Logfile.abort ()

       try :
          ret = self.extern.process ()

          if ret : msg = 'Program finished'
          else :   msg = 'Program finished with error'

       except KeyboardInterrupt :
          msg = 'Program aborted by Control C'; ret = False
       
       self.extern.finish ()      
       Logfile.showLabel (msg)
       return ret

#endclass

# -------------------------------------------------------------------------------------------------
#  Template for derived class of MainObj

class ExternMainObj (MainObj) :
    
    def __init__ (self) :

        MainObj.__init__ (self, self, version=xx, runTimeLog=None, errorLog=None)
        # own init operations

    def init (self) :
        # own operations
        return True

    def process (self) :
        #own operations
        return True

    def finish (self) :
        #own operations
        pass

#endclass
# -------------------------------------------------------------------------------------------------

def startTest (action, workingDir) :

    events   = 'NEAR-COAST-OF-NORTHERN-CHILE_2014-04-03T02-43-14'
    dir  = os.getcwd()
    dir  = Globals.setEventDir (os.path.join (dir,Globals.EVENTS, events))

    args = [sys.argv[0], action, '-f', dir]

    dir = Basic.changeDirectory (workingDir)      # create and/or change working director
    return args

# -------------------------------------------------------------------------------------------------
#   Module Test
# -------------------------------------------------------------------------------------------------

from time import sleep

VERSION_STRING = '<Version String>'

class TestObj (MainObj) :
    
    def __init__ (self, processFkt) :

        MainObj.__init__ (self, self, VERSION_STRING, 'TestProgram_run.log', 'TestProgram.log')
        self.processFkt = processFkt

    def init (self) :
        Logfile.add ('Init reached'); return True

    def process (self) :
        Logfile.add ('Process reached'); return self.processFkt()

    def finish (self) :
        Logfile.add ('Finish reached')

# -------------------------------------------------------------------------------------------------

def process1 () : 
    return True

def process2 () :

    Logfile.add ('Press Ctrl C')
    sleep (1000)
    return True

def testAll () :

    test1 = TestObj (process1)
    test1.run ()

    test2 = TestObj (process2)
    test2.run ()

    return True
