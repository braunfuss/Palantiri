from pyrocko.client import fdsn, catalog
import os
from pyrocko import util, io, trace, model, gf
import sys
sys.path.append('../tools/')
sys.path.append('../Common/')
import config
import Globals
import Basic
from optparse import OptionParser
from configparser import SafeConfigParser
from pyrocko import util, io, trace, cake
from config import Event, Trigger
from ConfigFile import ConfigObj, FilterCfg, OriginCfg
global options, args
from pyrocko.io import stationxml
from pyrocko.gf import ws, LocalEngine, Target, DCSource, RectangularSource
from pyrocko import util, model
from pyrocko.client import catalog
import numpy as num

def main(args):

    parser = OptionParser(usage="\npython %prog -t 2009-12-31T12:23:34 -d 5 -m SDS -k key -s all/acq")

    parser.add_option("-t","--time",      type="string", dest="time",      help="time")
    parser.add_option("-d","--duration",  type="string", dest="duration",  help="duration in min")
    parser.add_option("-m","--sdsfolder", type="string", dest="sdsfolder", help="sdsfolder")
    parser.add_option("-s","--station",   type="string", dest="stationversion", help="stationversion")
    parser.add_option("-f","--evpath",    type="string", dest="eventpath", help="eventpath")
    parser.add_option("-x","--dummy",     type="string", dest="station",   help="dummy")
    parser.add_option("-n","--dummy2",    type="string", dest="network",   help="dummy2")

    return parser.parse_args(args)


def globalConf():

    cDict = {}
    parser = SafeConfigParser()
    parser.read(os.path.join('../', 'global.conf'))

    for section_name in parser.sections():
        for name, value in parser.items(section_name):
            cDict[name] = value

    return cDict


params = globalConf()
options, args = main(sys.argv)
Basic.checkExistsDir(options.eventpath, isAbort=True)
Globals.setEventDir(options.eventpath)
C = config.Config(options.eventpath)
Origin = C.parseConfig('origin')
Conf = globalConf()
Config = C.parseConfig('config')

filter = FilterCfg(Config)

cfg = ConfigObj(dict=Config)
minDist = float(params['mindist'])
maxDist = float(params['maxdist'])
ev = Event(Origin['lat'], Origin['lon'], Origin['depth'], Origin['time'])
event = model.Event(lat=float(ev.lat), lon=float(ev.lon),
                    depth=float(ev.depth)*1000.,
                    time=util.str_to_time(ev.time))
tmin = util.str_to_time(ev.time)-40
tmax = util.str_to_time(ev.time)+40
global_cmt_catalog = catalog.GlobalCMT()
events = global_cmt_catalog.get_events(
    time_range=(tmin, tmax),
    latmin=float(ev.lat)-1.,
    latmax=float(ev.lat)+1,
    lonmin=float(ev.lon)-1,
    lonmax=float(ev.lon)+1)
event_cat = events[0]


source = gf.DCSource(lat=event_cat.lat, lon=event_cat.lon,
                     strike=event_cat.moment_tensor.strike1,
                     rake=event_cat.moment_tensor.rake1,
                     dip=event_cat.moment_tensor.dip1,
                     magnitude=event.magnitude)
newFreq = float(filter.newFrequency())
options.time = Origin['time']
options.duration = int(Conf['duration'])
sdspath = os.path.join(options.eventpath, 'data')
try:
    os.mkdir(sdspath)
except Exception:
    pass
model.dump_events([event], sdspath+'event.pf')
if float(params['duration']) == 0:
    tmin = util.str_to_time(ev.time)+float(params['tmin'])
    tmax = util.str_to_time(ev.time)+float(params['tmax'])
else:
    tmin = util.str_to_time(ev.time)
    tmax = util.str_to_time(ev.time) + float(params['duration'])

def get_stations(site, lat, lon, rmin, rmax, tmin, tmax,
                 channel_pattern='BH*'):
    extra = {}
    if site == 'iris':
        extra.update(matchtimeseries=True)

    sx = fdsn.station(
            site=site, latitude=lat, longitude=lon,
            minradius=rmin, maxradius=rmax,
            startbefore=tmin, endafter=tmax, channel=channel_pattern,
            format='text', level='channel', includerestricted=False)

    return sx.get_pyrocko_stations()


site = 'geofon'
minDist = float(params['mindist'])
maxDist = float(params['maxdist'])
diffDist = (maxDist - minDist)/9.
displacement_geofon = []
stations_disp_geofon = []
stations_real_geofon = []
gaps = []
trs_projected_geofon = []
trs_projected_displacement_geofon = []

quantity_to_unit = {
    'displacement': 'M',
    'velocity': 'M/S',
    'acceleration': 'M/S**2'}

quantity = cfg.quantity()


try:
    for l in range(0, 1):
        gaps = []
        stations_geofon = get_stations(site, event.lat, event.lon, minDist,
                                       maxDist, tmin, tmax, 'BH*')

        nstations_geofon = [s for s in stations_geofon]

        selection_geofon = fdsn.make_data_selection(nstations_geofon, tmin,
                                                    tmax)
        request_waveform_geofon = fdsn.dataselect(site=site,
                                                  selection=selection_geofon)

        with open(os.path.join(sdspath,
                               'traces_geofon_part%s.mseed' % l), 'wb') as file:
            file.write(request_waveform_geofon.read())
        print('traces written')
        traces_geofon = io.load(os.path.join(sdspath,
                                             'traces_geofon_part%s.mseed' % l))

        for tr in traces_geofon:
            for st in stations_geofon:
                for channel in st.channels:
                    if tr.station == st.station and tr.location == st.location and channel.name == tr.channel and tr.location == st.location and tr.network == st.network:
                        stations_real_geofon.append(st)
                        gaps.append(st.station)
        remove = [x for x in gaps if gaps.count(x) > 3]
        for re in remove:
            stations_real_geofon.remove(re)

        request_response = fdsn.station(
            site=site, selection=selection_geofon, level='response')
        request_response.dump(filename=os.path.join(sdspath,
                              'responses_geofon_part%s.yml' % l))
        request_response.dump_xml(filename=os.path.join(sdspath,
                                  'responses_geofon_part%s.xml' % l))
        sx = stationxml.load_xml(filename=os.path.join(sdspath,
                                 'responses_geofon_part%s.xml' % l))
        pyrocko_stations = sx.get_pyrocko_stations()
        event_origin = gf.Source(
                                lat=event.lat,
                                lon=event.lon)

        traces_geofon = io.load(os.path.join(sdspath,
                                             'traces_geofon_part%s.mseed' % l))

        projections = []
        for station in stations_real_geofon:

                    backazimuth = source.azibazi_to(station)[1]
                    projections.extend(station.guess_projections_to_rtu(
                                        out_channels=('R', 'T', 'Z'),
                                        backazimuth=backazimuth))


        for matrix, in_channels, out_channels in projections:
            deps = trace.project_dependencies(
            matrix, in_channels, out_channels)

        trs_projected_geofon.extend(
                trace.project(
                traces_geofon, matrix,
                in_channels, out_channels))

        disp_rot = []
        for tr in traces_geofon:
            for station in stations_real_geofon:
                if tr.station == station.station and tr.location == station.location:
                    try:
                        polezero_response = request_response.get_pyrocko_response(
                        nslc=tr.nslc_id,
                        timespan=(tr.tmin, tr.tmax),
                        fake_input_units=quantity_to_unit[quantity])

                        restituted = tr.transfer(
                        tfade=2.,
                        freqlimits=(0.01, 0.1, 1., 2.),
                        transfer_function=polezero_response,
                        invert=True)
                        if quantity == 'velocity':
                            tr.ydata = num.diff(tr.ydata)
                        displacement_geofon.append(restituted)
                        disp_rot.append(restituted)
                        stations_disp_geofon.append(station)
                    except:
                        pass
        model.dump_stations(stations_disp_geofon,
                            os.path.join(sdspath,
                                         'stations_geofon_part%s.txt' % l))
        trs_projected_displacement_geofon.extend(
                trace.project(
                disp_rot, matrix,
                in_channels, out_channels))
except:
    for l in range(0,9):
        try:
            gaps = []
            maxDist = minDist+diffDist
            stations_geofon = get_stations(site, event.lat,event.lon,minDist, maxDist,tmin,tmax, 'BH*')

            nstations_geofon = [s for s in stations_geofon]

            selection_geofon = fdsn.make_data_selection(nstations_geofon, tmin, tmax)
            request_waveform_geofon = fdsn.dataselect(site=site, selection=selection_geofon)

            with open(os.path.join(sdspath,'traces_geofon_part%s.mseed' %l), 'wb') as file:
                file.write(request_waveform_geofon.read())
            print('traces written')
            traces_geofon = io.load(os.path.join(sdspath,'traces_geofon_part%s.mseed' %l))

            for tr in traces_geofon:
                for st in stations_geofon:
                    for channel in st.channels:
                        if tr.station == st.station and tr.location == st.location and channel.name == tr.channel and tr.location == st.location and tr.network == st.network:
                            stations_real_geofon.append(st)
                            gaps.append(st.station)
            remove = [x for x in gaps if gaps.count(x) > 3]
            for re in remove:
                stations_real_geofon.remove(re)
            request_response = fdsn.station(
                site=site, selection=selection_geofon, level='response')
            request_response.dump(filename=os.path.join(sdspath,'responses_geofon_part%s.yml'%l))
            request_response.dump_xml(filename=os.path.join(sdspath,'responses_geofon_part%s.xml'%l))
            sx = stationxml.load_xml(filename=os.path.join(sdspath,'responses_geofon_part%s.xml'%l))
            pyrocko_stations = sx.get_pyrocko_stations()
            event_origin = gf.Source(
            lat=event.lat,
            lon=event.lon)

            traces_geofon = io.load(os.path.join(sdspath,'traces_geofon_part%s.mseed' %l))

            projections = []
            for station in stations_real_geofon:
                        backazimuth = source.azibazi_to(station)[1]
                        projections.extend(station.guess_projections_to_rtu(
                        out_channels=('R', 'T', 'Z'),
                        backazimuth=backazimuth))


            for matrix, in_channels, out_channels in projections:
                deps = trace.project_dependencies(
                matrix, in_channels, out_channels)

            trs_projected_geofon.extend(
                    trace.project(
                    traces_geofon, matrix,
                    in_channels, out_channels))
            disp_rot = []
            for tr in traces_geofon:
                for station in stations_real_geofon:
                    if tr.station == station.station and tr.location == station.location:
                        try:
                            polezero_response = request_response.get_pyrocko_response(
                            nslc=tr.nslc_id,
                            timespan=(tr.tmin, tr.tmax),
                            fake_input_units=quantity_to_unit[quantity])

                            restituted = tr.transfer(
                            tfade=2.,
                            freqlimits=(0.01, 0.1, 1., 2.),
                            transfer_function=polezero_response,
                            invert=True)
                            if quantity == 'velocity':
                                tr.ydata = num.diff(tr.ydata)
                            disp_rot.append(restituted)

                            displacement_geofon.append(restituted)
                            stations_disp_geofon.append(station)
                        except:
                            pass
        except:
            pass
        model.dump_stations(stations_disp_geofon, os.path.join(sdspath,'stations_geofon_part%s.txt' %l))

        minDist = minDist+diffDist
        try:
            trs_projected_displacement_geofon.extend(
                    trace.project(
                    disp_rot, matrix,
                    in_channels, out_channels))
        except:
            pass

io.save(displacement_geofon, os.path.join(sdspath,'traces_restituted_geofon.mseed'))
model.dump_stations(stations_disp_geofon, os.path.join(sdspath,'stations_disp_geofon.txt'))
model.dump_stations(stations_disp_geofon, os.path.join(sdspath,'stations_geofon.txt'))
io.save(trs_projected_displacement_geofon, os.path.join(sdspath,'traces_restituted_rotated_geofon.mseed'))
io.save(trs_projected_geofon, os.path.join(sdspath,'traces_rotated_geofon.mseed'))



stations_sites = []
stations_real_sites = []
stations_disp_sites = []
displacement_sites = []
traces_sites = []
trs_projected = []
trs_projected_displacement = []
minDist = float(params['mindist'])
maxDist = float(params['maxdist'])

sites = ['iris','orfeus', 'resif', 'usp', 'bgr', 'ingv', 'geonet', 'ethz', 'ncedc', 'knmi', 'isc', 'ipgp', 'koeri']
for site in sites:
    try:
        stations_site = get_stations(site, event.lat,event.lon,minDist, maxDist,tmin,tmax, 'BH*')
        if not stations_sites:
            stations_sites = stations_site
        else:
            stations_sites = stations_sites + stations_site

        nstations_site = [s for s in stations_site]

        selection_site = fdsn.make_data_selection(nstations_site, tmin, tmax)
        request_waveform_site = fdsn.dataselect(site=site, selection=selection_site)
        with open(os.path.join(sdspath,'traces_%s.mseed' %site), 'wb') as file:
            file.write(request_waveform_site.read())
        print('traces written')
        traces_site = io.load(os.path.join(sdspath,'traces_%s.mseed'%site))
        stations_real_site = []

        if not traces_sites:
            traces_sites = traces_site
        else:
            traces_sites = traces_sites+traces_site
        gaps= []

        for tr in traces_site:
            for st in stations_site:
                for channel in st.channels:
                    if tr.station == st.station and tr.location == st.location and channel.name == tr.channel and tr.location == st.location and tr.network == st.network:
                        stations_real_site.append(st)
                        stations_real_sites.append(st)
                        gaps.append(st.station)
        remove = [x for x in gaps if gaps.count(x) > 3]
        for re in remove:
                stations_real_site.remove(re)

        request_response = fdsn.station(
            site=site, selection=selection_site, level='response')
        request_response.dump(filename=os.path.join(sdspath,'responses_%s.yml'%site))
        request_response.dump_xml(filename=os.path.join(sdspath,'responses_%s.xml'%site))
        sx = stationxml.load_xml(filename=os.path.join(sdspath,'responses_%s.xml'%site))
        pyrocko_stations = sx.get_pyrocko_stations()
        event_origin = gf.Source(
        lat=event.lat,
        lon=event.lon)

        traces_site = io.load(os.path.join(sdspath,'traces_%s.mseed'%site))

        projections = []
        for station in stations_real_site:
                    backazimuth = source.azibazi_to(station)[1]
                    projections.extend(station.guess_projections_to_rtu(
                    out_channels=('R', 'T', 'Z'),
                    backazimuth=backazimuth))

        for matrix, in_channels, out_channels in projections:
            deps = trace.project_dependencies(
            matrix, in_channels, out_channels)

        trs_projected.extend(
                trace.project(
                traces_site, matrix,
                in_channels, out_channels))

        displacement_site = []
        stations_disp_site = []
        disp_rot = []
        for tr in traces_site:
            for station in stations_real_site:
                if tr.station == station.station and tr.location == station.location:
                    try:
                        polezero_response = request_response.get_pyrocko_response(
                        nslc=tr.nslc_id,
                        timespan=(tr.tmin, tr.tmax),
                        fake_input_units=quantity_to_unit[quantity])

                        restituted = tr.transfer(
                        tfade=2.,
                        freqlimits=(0.01, 0.1, 1., 2.),
                        transfer_function=polezero_response,
                        invert=True)
                        if quantity == 'velocity':
                            tr.ydata = num.diff(tr.ydata)

                        displacement_site.append(restituted)
                        displacement_sites.append(restituted)
                        disp_rot.append(restituted)
                        stations_disp_site.append(station)
                        stations_disp_sites.append(station)
                    except:
                        pass
        trs_projected_displacement.extend(
                trace.project(
                disp_rot, matrix,
                in_channels, out_channels))
        io.save(displacement_site, os.path.join(sdspath,'traces_restituted_%s.mseed'%site))
        model.dump_stations(stations_disp_site, os.path.join(sdspath,'stations_disp_%s.txt'%site))
        model.dump_stations(stations_disp_site, os.path.join(sdspath,'stations_%s.txt'%site))

    except:
        pass

stations_all  = stations_disp_sites+stations_disp_geofon
for stg in stations_real_geofon:
    for sti in stations_real_sites:
        if sti.station == stg.station and sti.location == stg.location and sti.network == stg.network:
            try:
                stations_all.remove(sti)
            except:
                pass
        else:
            pass
try:
    traces_all = traces_sites+traces_geofon
    traces_all_rot = trs_projected_geofon+trs_projected
except:
    traces_all = traces_sites
    traces_all_rot = trs_projected

for tr in traces_all:
    try:
        tr.downsample_to(newFreq)
    except:
        pass
io.save(traces_all, os.path.join(sdspath,'traces.mseed'))
model.dump_stations(stations_all, os.path.join(sdspath,'stations.txt'))
io.save(trs_projected_displacement, os.path.join(sdspath,'traces_restituted_rotated_sites.mseed'))
io.save(trs_projected, os.path.join(sdspath,'traces_rotated_sites.mseed'))
io.save(traces_all_rot, os.path.join(sdspath,'traces_rotated.mseed'))


try:
    stations_all_disp = stations_disp_sites+stations_disp_geofon
    for stg in stations_disp_geofon:
        for sti in stations_disp_sites:
            if sti.station == stg.station and sti.location == stg.location:
                try:
                    stations_all_disp.remove(sti)
                except:
                    pass
            else:
                pass
    traces_all_disp = displacement_sites+displacement_geofon
    traces_all_rot_disp = trs_projected_displacement_geofon+trs_projected_displacement
    for tr in traces_all_disp:
        tr.downsample_to(newFreq)
    for tr in traces_all_rot_disp:
        tr.downsample_to(newFreq)

    for st in stations_all_disp:
        for channel in st.channels:
            if channel.name == 'BHE':
                channel.name = 'R'
            if channel.name == 'BHN':
                channel.name = 'T'
            if channel.name == 'BHZ':
                channel.name = 'Z'
    io.save(traces_all_rot_disp, os.path.join(sdspath,'traces_restituted_rotated.mseed'))
    io.save(traces_all_disp, os.path.join(sdspath,'traces_restituted.mseed'))
    model.dump_stations(stations_all_disp, os.path.join(sdspath,'stations_disp.txt'))
except:
    stations_all_disp = stations_disp_sites
    traces_all_disp = displacement_sites
    traces_all_rot_disp = trs_projected_displacement
    for tr in traces_all_disp:
        tr.downsample_to(newFreq)
    for tr in traces_all_rot_disp:
        tr.downsample_to(newFreq)
    for st in stations_all_disp:
        for channel in st.channels:
            if channel.name == 'BHE':
                channel.name = 'R'
            if channel.name == 'BHN':
                channel.name = 'T'
            if channel.name == 'BHZ':
                channel.name = 'Z'
    io.save(traces_all_rot_disp, os.path.join(sdspath,'traces_restituted_rotated.mseed'))
    io.save(traces_all_disp, os.path.join(sdspath,'traces_restituted.mseed'))
    model.dump_stations(stations_all_disp, os.path.join(sdspath,'stations_disp.txt'))
