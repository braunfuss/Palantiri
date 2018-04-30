
import os
import sys
import imp

import Basic, Globals                        # Import from Common

import getStationList, getStationWaveformData, ev_meta_mt4, WaveformProgs

# -------------------------------------------------------------------------------------------------

Basic.onlyForInternalUse ()

action = 'getstations'
#action = 'getdata'
action = 'getmeta'

events = 'NEAR-COAST-OF-NORTHERN-CHILE_2014-04-03T02-43-14'
dir    = Basic.changeDirectory ('..')
#dir   = Globals.setEventDir   (os.path.join (dir, Globals.EVENTS, events))

if len (sys.argv) == 1 :
   if action == 'getstations' : sys.argv = [sys.argv[0], action, events]
   else :
      #sys.argv = [sys.argv[0], action, events, 'geofon']
       sys.argv = [sys.argv[0], action, events, 'iris']
      #sys.argv = [sys.argv[0], action, events, 'AK']
     
   config = imp.load_source ('Config',os.path.join (os.getcwd(),'tools','config.py'))
   cmd    = WaveformProgs.buildCommand (config)

   oldArgs  = sys.argv
   sys.argv = cmd.split()
#endif

dir = Basic.changeDirectory ('tools')      # create and/or change working directory

if   action == 'getstations' : getStationList.MainProc()
elif action == 'getdata' :     getStationWaveformData.MainProc()
elif action == 'getmeta' :     ev_meta_mt4.MainProc()
else :                         assert False

#print '*** End of test ***'
