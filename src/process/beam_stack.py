from pyrocko.model import Station, dump_stations
from pyrocko.guts import Object, Float, String, Bool, Dict
from pyrocko import orthodrome as ortho
from pyrocko import util, io, trace, cake
import logging
import numpy as num
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from collections import defaultdict
from matplotlib import cm

logger = logging.getLogger('beam-forming')


r_earth = 6371000.785
torad = num.pi/180.
onedeg = r_earth*torad


def to_cartesian(items, reflatlon):
    res = defaultdict()
    for i, item in enumerate(items):

        y, x = ortho.latlon_to_ne(reflatlon, item)
        depth = item.depth
        elevation = item.elevation
        dz = elevation - depth
        lat = item.lat/180.*num.pi
        z = r_earth+dz*num.sin(lat)
        res[item.nsl()[:2]] = (x, y, z)
    return res


def cmp(a, b):
    return (a > b) - (a < b)


class BeamForming(Object):
    station_c = Station.T(optional=True)
    bazi = Float.T()
    slow = Float.T()
    diff_dt_treat = String.T(help='how to handle differing sampling rates:'
                             ' oversample(default) or downsample')
    normalize_std = Bool.T()
    post_normalize = Bool.T()
    t_shifts = Dict.T(String.T(), Float.T())

    def __init__(self, stations, traces, normalize=True, post_normalize=False,
                 diff_dt_treat='oversample'):
        self.stations = stations
        self.c_lat_lon_z = self.center_lat_lon(stations)
        self.traces = traces
        self.diff_dt_treat = diff_dt_treat
        self.normalize_std = normalize
        self.post_normalize = post_normalize
        self.station_c = None
        self.diff_dt_treat = diff_dt_treat

    def process(self, event, timing, bazi=None, slow=None,  restitute=False,
                *args, **kwargs):
        '''
      :param timing: CakeTiming. Uses the definition without the offset.
      :param fn_dump_center: filename to where center stations shall be dumped
      :param fn_beam: filename of beam trace
      :param model: earthmodel to use(optional)
      :param earthmodel to use(optional)
      :param network: network code(optional)
      :param station: station code(optional)
        '''
        logger.debug('start beam forming')
        stations = self.stations
        network_code = kwargs.get('responses', None)
        network_code = kwargs.get('network', '')
        station_code = kwargs.get('station', 'STK')
        c_station_id = (network_code, station_code)
        t_shifts = []
        lat_c, lon_c, z_c = self.c_lat_lon_z

        self.station_c = Station(lat=float(lat_c),
                                 lon=float(lon_c),
                                 elevation=float(z_c),
                                 depth=0.,
                                 name='Array Center',
                                 network=c_station_id[0],
                                 station=c_station_id[1][:5])
        fn_dump_center = kwargs.get('fn_dump_center', 'array_center.pf')
        fn_beam = kwargs.get('fn_beam', 'beam.mseed')
        if event:
            mod = cake.load_model(crust2_profile=(event.lat, event.lon))
            dist = ortho.distance_accurate50m(event, self.station_c)
            ray = timing.t(mod, (event.depth, dist), get_ray=True)

            if ray is None:
                logger.error('None of defined phases available at beam \
                              station:\n %s' % self.station_c)
                return
            else:
                b = ortho.azimuth(self.station_c, event)
                if b>=0.:
                    self.bazi = b
                elif b<0.:
                    self.bazi = 360.+b
                self.slow = ray.p/(cake.r2d*cake.d2m)
        else:
            self.bazi = bazi
            self.slow = slow

        logger.info('stacking %s with slowness %1.4f s/km at back azimut %1.1f'
                    'degrees' % ('.'.join(c_station_id),
                                 self.slow*cake.km, self.bazi))

        lat0 = num.array([lat_c]*len(stations))
        lon0 = num.array([lon_c]*len(stations))
        lats = num.array([s.lat for s in stations])
        lons = num.array([s.lon for s in stations])
        ns, es = ortho.latlon_to_ne_numpy(lat0, lon0, lats, lons)
        theta = num.float(self.bazi*num.pi/180.)
        R = num.array([[num.cos(theta), -num.sin(theta)],
                       [num.sin(theta), num.cos(theta)]])
        distances = R.dot(num.vstack((es, ns)))[1]
        channels = set()
        self.stacked = {}
        num_stacked = {}
        self.t_shifts = {}
        self.shifted_traces = []
        taperer = trace.CosFader(xfrac=0.05)
        if self.diff_dt_treat == 'downsample':
            self.traces.sort(key=lambda x: x.deltat)
        elif self.diff_dt_treat == 'oversample':
            dts = [t.deltat for t in self.traces]
            for tr in self.traces:
                tr.resample(min(dts))

        for tr in self.traces:
            if tr.nslc_id[:2] == c_station_id:
                continue
            tr = tr.copy(data=True)
            tr.ydata = tr.ydata.astype(num.float64) - tr.ydata.mean(dtype=num.float64)
            tr.taper(taperer)
            try:
                stack_trace = self.stacked[tr.channel]
                num_stacked[tr.channel] += 1
            except KeyError:
                stack_trace = tr.copy(data=True)
                stack_trace.set_ydata(num.zeros(
                    len(stack_trace.get_ydata())))

                stack_trace.set_codes(network=c_station_id[0],
                                      station=c_station_id[1],
                                      location='',
                                      channel=tr.channel)

                self.stacked[tr.channel] = stack_trace
                channels.add(tr.channel)
                num_stacked[tr.channel] = 1

            nslc_id = tr.nslc_id

            try:
                stats = list(filter(lambda x: util.match_nslc(
                    '%s.%s.%s.*' % x.nsl(), nslc_id), stations))
                stat = stats[0]
            except IndexError:
                break

            i = stations.index(stat)
            d = distances[i]
            t_shift = d*self.slow
            t_shifts.append(t_shift)
            tr.shift(t_shift)
            self.t_shifts[tr.nslc_id[:2]] = t_shift
            if self.normalize_std:
                tr.ydata = tr.ydata/tr.ydata.std()

            if num.abs(tr.deltat-stack_trace.deltat) > 0.000001:
                if self.diff_dt_treat == 'downsample':
                    stack_trace.downsample_to(tr.deltat)
                elif self.diff_dt_treat == 'upsample':
                    raise Exception('something went wrong with the upsampling,\
                                     previously')
            stack_trace.add(tr)
            self.shifted_traces.append(tr)

        if self.post_normalize:
            for ch, tr in self.stacked.items():
                tr.set_ydata(tr.get_ydata()/num_stacked[ch])

        self.save_station(fn_dump_center)
        self.checked_nslc([stack_trace])
        self.save(stack_trace, fn_beam)
        return self.shifted_traces, stack_trace, t_shifts

    def checked_nslc(self, trs):
        for tr in trs:
            oldids = tr.nslc_id
            n, s, l, c = oldids
            tr.set_codes(network=n[:2], station=s[:5], location=l[:2],
                         channel=c[:3])
            newids = tr.nslc_id
            if cmp(oldids, newids) != 0:
                logger.warn('nslc id truncated: %s to %s' %
                            ('.'.join(oldids), '.'.join(newids)))

    def snuffle(self):
        '''Scrutinize the shifted traces.'''
        from pyrocko import snuffler
        snuffler.snuffle(self.shifted_traces)

    def center_lat_lon(self, stations):
        '''Calculate a mean geographical centre of the array
        using spherical earth'''

        lats = num.zeros(len(stations))
        lons = num.zeros(len(stations))
        elevations = num.zeros(len(stations))
        depths = num.zeros(len(stations))
        for i, s in enumerate(stations):
            lats[i] = s.lat*torad
            lons[i] = s.lon*torad
            depths[i] = s.depth
            elevations[i] = s.elevation

        z = num.mean(elevations-depths)
        return(lats.mean()*180/num.pi, lons.mean()*180/num.pi, z)

    def plot(self, fn='beam_shifts.png'):
        stations = self.stations
        stations.append(self.station_c)
        res = to_cartesian(stations, self.station_c)
        center_xyz = res[self.station_c.nsl()[:2]]
        x = num.zeros(len(res))
        y = num.zeros(len(res))
        z = num.zeros(len(res))
        sizes = num.zeros(len(res))
        stat_labels = []
        i = 0
        for nsl, xyz in res.items():
            x[i] = xyz[0]
            y[i] = xyz[1]
            z[i] = xyz[2]

            try:
                sizes[i] = self.t_shifts[nsl[:2]]
                stat_labels.append('%s' % ('.'.join(nsl)))
            except AttributeError:
                self.fail('Run the snuffling first')
            except KeyError:
                stat_labels.append('%s' % ('.'.join(nsl)))
                continue
            finally:
                i += 1

        x /= 1000.
        y /= 1000.
        z /= 1000.
        xmax = x.max()
        xmin = x.min()
        ymax = y.max()
        ymin = y.min()

        fig = plt.figure()
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax = fig.add_axes([0.10, 0.1, 0.70, 0.7])
        ax.set_aspect('equal')
        cmap = cm.get_cmap('bwr')
        ax.scatter(x, y, c=sizes, s=200, cmap=cmap,
                   vmax=num.max(sizes), vmin=-num.max(sizes))
        for i, lab in enumerate(stat_labels):
            ax.text(x[i], y[i], lab, size=14)

        x = x[num.where(sizes == 0.)]
        y = y[num.where(sizes == 0.)]
        ax.scatter(x, y, c='black', alpha=0.4, s=200)

        ax.arrow(center_xyz[0]/1000.,
                 center_xyz[1]/1000.,
                 -num.sin(self.bazi/180.*num.pi),
                 -num.cos(self.bazi/180.*num.pi),
                 head_width=0.2,
                 head_length=0.2)
        ax.set_ylabel("N-S [km]")
        ax.set_xlabel("E-W [km]")
        ColorbarBase(cax, cmap=cmap,
                     norm=Normalize(vmin=sizes.min(), vmax=sizes.max()))
        logger.debug('finished plotting %s' % fn)
        fig.savefig(fn)

    def save(self, traces, fn='beam.pf'):
        io.save(traces, fn)

    def save_station(self, fn):
        dump_stations([self.station_c], fn)
