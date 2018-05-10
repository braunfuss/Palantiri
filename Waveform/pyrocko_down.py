from pyrocko.client import fdsn, catalog
import os
from pyrocko import util, io, trace, model, gf
import sys
sys.path.append ('../tools/')
sys.path.append ('../Common/')
import config
import Globals                     # Own global data
import Basic                       # Own module with basic functions
import Logfile                     # Implements logfile
import ConfigFile                  # Semantic of config file entries
import Debug
import DataTypes                   # Common data types
import KeyFile
import DataDir
from optparse import OptionParser
from   ConfigParser import SafeConfigParser
from pyrocko import util, io, trace, cake
from   config import Event, Trigger
from    ConfigFile import ConfigObj, FilterCfg, OriginCfg
global options,args



def main (args):

    parser = OptionParser (usage="\npython %prog -t 2009-12-31T12:23:34 -d 5 -m SDS -k key -s all/acq")

    parser.add_option ("-t","--time",      type="string", dest="time",      help="time")
    parser.add_option ("-d","--duration",  type="string", dest="duration",  help="duration in min")
    parser.add_option ("-m","--sdsfolder", type="string", dest="sdsfolder", help="sdsfolder")
    parser.add_option ("-s","--station",   type="string", dest="stationversion", help="stationversion")
    parser.add_option ("-f","--evpath",    type="string", dest="eventpath", help="eventpath")
    parser.add_option ("-x","--dummy",     type="string", dest="station",   help="dummy")    #hs : client flag
    parser.add_option ("-n","--dummy2",    type="string", dest="network",   help="dummy2")   #hs : single network

    return parser.parse_args(args)

def globalConf():

    cDict  = {}
    parser = SafeConfigParser()
    parser.read (os.path.join ('..', 'global.conf'))

    for section_name in parser.sections():
        for name, value in parser.items(section_name) : cDict[name] = value

    return cDict

options,args = main (sys.argv)
Basic.checkExistsDir (options.eventpath, isAbort=True)
Globals.setEventDir  (options.eventpath)
C      = config.Config (options.eventpath)
Origin = C.parseConfig ('origin')
Conf   = globalConf()
Config = C.parseConfig ('config')

filter = FilterCfg (Config)

cfg     = ConfigObj (dict=Config)
minDist, maxDist = cfg.FloatRange ('mindist', 'maxdist')

ev = Event (Origin['lat'],Origin['lon'],Origin['depth'],Origin['time'] )
event = model.Event(lat=float(ev.lat), lon=float(ev.lon), depth=float(ev.depth)*1000., time=util.str_to_time(ev.time))
newFreq                = float (filter.newFrequency())
options.time           = Origin ['time']
options.duration       = int (Conf['duration'])
#options.sdsfolder      = os.path.join (options.eventpath,'data')
sdspath = os.path.join(options.eventpath,'data')
tmin = util.str_to_time(ev.time)-600.
tmax = util.str_to_time(ev.time)+1800.
site = 'iris'

def get_stations(site, lat, lon, rmin, rmax, tmin, tmax, channel_pattern='BH*'):
    from pyrocko.fdsn import ws
    extra = {}
    if site == 'iris':
        extra.update(matchtimeseries=True)

    sx = fdsn.station(
		site=site, latitude=lat, longitude=lon,
		minradius=rmin, maxradius=rmax,
		startbefore=tmin, endafter=tmax, channel=channel_pattern,
		format='text', level='channel', includerestricted=False)

    return sx.get_pyrocko_stations()


stations = get_stations(site, event.lat,event.lon,minDist, maxDist,tmin,tmax, 'BHZ')

#model.dump_stations(stations, os.path.join (sdspath,'stations.txt'))
# setup a waveform data request

nstations = [s for s in stations]


selection = fdsn.make_data_selection(nstations, tmin, tmax)
request_waveform = fdsn.dataselect(site=site, selection=selection)

# write the incoming data stream to 'traces.mseed'
with open(os.path.join (sdspath,'traces.mseed'), 'wb') as file:
    file.write(request_waveform.read())
print('traces written')
# request meta data
traces = io.load(os.path.join (sdspath,'traces.mseed'))

stations_real = []
gaps= []
for tr in traces:
    for st in stations:
        if tr.station == st.station and tr.location == st.location:
                stations_real.append(st)
                gaps.append(st.station)
remove =[x for x in gaps if gaps.count(x) > 1]
for re in remove:
    for st in stations_real:
        if st.station == re:
            stations_real.remove(st)
model.dump_stations(stations_real, os.path.join (sdspath,'stations.txt'))

request_response = fdsn.station(
    site=site, selection=selection, level='response')
from pyrocko.io import stationxml
# save the response in YAML and StationXML format
request_response.dump(filename=os.path.join (sdspath,'responses.yml'))
request_response.dump_xml (filename=os.path.join (sdspath,'responses.xml'))
sx = stationxml.load_xml(filename=os.path.join (sdspath,'responses.xml'))
pyrocko_stations = sx.get_pyrocko_stations()
model.dump_stations(stations_real, os.path.join (sdspath,'stations2.txt'))

# Loop through retrieved waveforms and request meta information
# for each trace
event_origin = gf.Source(
lat=event.lat,
lon=event.lon)

traces = io.load(os.path.join (sdspath,'traces.mseed'))


displacement = []
for tr in traces:
    try:
        polezero_response = request_response.get_pyrocko_response(
        nslc=tr.nslc_id,
        timespan=(tr.tmin, tr.tmax),
        fake_input_units='M')
        # *fake_input_units*: required for consistent responses throughout entire
        # data set

        # deconvolve transfer function
        restituted = tr.transfer(
        tfade=2.,
        freqlimits=(0.01, 0.1, 1., 2.),
        transfer_function=polezero_response,
        invert=True)

        displacement.append(restituted)
    except:
        pass


io.save(displacement, os.path.join (sdspath,'traces_restituted.mseed'))
