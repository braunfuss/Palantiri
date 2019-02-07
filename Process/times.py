import sys
sys.path.append ('../Common/')

from obspy.core.utcdatetime import UTCDateTime

import Logfile

def calculateTimeWindows(mint,maxt,Config,Origin, switch):

    tw = {}
    st = str(Origin.time)[:-1]

    if switch == 0:
        winlen = float(Config['winlen'])
    if switch == 1:
        winlen = float(Config['winlen_f2'])

    tw['start'] = UTCDateTime(UTCDateTime(st)+(mint-float(Config['forerun'])))
    tw['end'] = tw['start']+float(Config['duration'])
    #tw['end']   = UTCDateTime(UTCDateTime(st)+(maxt+float(Config['duration'])+winlen))
    #timespan = UTCDateTime(UTCDateTime(st)+(maxt+float(Config['duration'])+winlen)) - UTCDateTime(UTCDateTime(st)+(mint-int(Config['forerun'])))
    timespan = tw['end']-tw['start']
    Logfile.red ('ORIGIN TIME %s ' % UTCDateTime(st))
    Logfile.red ('TIME WINDOW: %s - %s ' % (tw['start'], tw['end']) )
    Logfile.red ('TIME SPAN: %s Minutes ' % (timespan/60))

    return tw
