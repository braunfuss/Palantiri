
#import logging
import sys
import platform

WINDOWS = (platform.system() == 'Windows')

# add local directories to import path                  
sys.path.append ('../Common/')

from obspy.core.utcdatetime import UTCDateTime

#logger = logging.getLogger(sys.argv[0])
import Logfile

def calculateTimeWindows(mint,maxt,Config,Origin):

    tw = {}
    st = str(Origin.time)[:-1]
    
    #tw['start'] = UTCDateTime(UTCDateTime(st)+(mint-int(Config['forerun'])-int(Config['security'])))
    #tw['end'] = UTCDateTime(UTCDateTime(st)+(maxt+int(Config['duration'])+int(Config['winlen'])+int(Config['security'])))
    
    tw['start'] = UTCDateTime(UTCDateTime(st)+(mint-int(Config['forerun'])))
    tw['end']   = UTCDateTime(UTCDateTime(st)+(maxt+int(Config['duration'])+int(Config['winlen'])))

    timespan = UTCDateTime(UTCDateTime(st)+(maxt+int(Config['duration'])+int(Config['winlen']))) - UTCDateTime(UTCDateTime(st)+(mint-int(Config['forerun'])))
    
    #logger.info ('\033[31m ORIGIN TIME %s \033[0m'% UTCDateTime(st))
    #logger.info ('\033[31m TIME WINDOW: %s - %s \033[0m' % (tw['start'], tw['end']) )
    #logger.info( '\033[31m TIME SPAN: %s \033[0m' % (timespan/60))

    Logfile.red ('ORIGIN TIME %s ' % UTCDateTime(st))
    Logfile.red ('TIME WINDOW: %s - %s ' % (tw['start'], tw['end']) )
    Logfile.red ('TIME SPAN: %s ' % (timespan/60))

    return tw



