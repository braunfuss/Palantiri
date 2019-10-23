from affine import Affine
from pyproj import Proj, transform
from mpl_toolkits.basemap import Basemap
import numpy as num
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from matplotlib.pyplot import cm
from matplotlib.widgets import Slider
from pylab import plot, show, figure, scatter, axes, draw
import csv
from obspy.imaging.beachball import beach
import matplotlib.colors as colors
import matplotlib.tri as tri
from pyrocko import trace, io, model, orthodrome
from matplotlib.colors import LinearSegmentedColormap
from palantiri.common import ConfigFile
from palantiri.common.ConfigFile import ConfigObj, FilterCfg, OriginCfg, SynthCfg
from palantiri.common import Globals
from palantiri.tools import config
from palantiri.process.sembCalc import toAzimuth
global evpath
from pyrocko import util

w=25480390.0


def draw_beach(ax, scale, map, event):
    desired=[3,4]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
    desired=[7,8,9]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
    x, y = map(event_cor[1][0],event_cor[0][0])
    ax = plt.gca()

    np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
    beach1 = beach(np1, xy=(x, y), width=0.09)
    ax.add_collection(beach1)


def draw_sources(ax, syn_in, map, scale):
    from pyrocko.gf import DCSource, RectangularSource, MTSource

    sources = []

    if syn_in.source() == 'RectangularSource':
        sources.append(RectangularSource(
            lat=float(syn_in.lat_0()),
            lon=float(syn_in.lon_0()),
            east_shift=float(syn_in.east_shift_0())*1000.,
            north_shift=float(syn_in.north_shift_0())*1000.,
            depth=syn_in.depth_syn_0()*1000.,
            strike=syn_in.strike_0(),
            dip=syn_in.dip_0(),
            rake=syn_in.rake_0(),
            width=syn_in.width_0()*1000.,
            length=syn_in.length_0()*1000.,
            nucleation_x=syn_in.nucleation_x_0(),
            slip=syn_in.slip_0(),
            nucleation_y=syn_in.nucleation_y_0(),
            velocity=syn_in.velocity_0(),
            anchor=syn_in.anchor(),
            time=util.str_to_time(syn_in.time_0())))
        for source in sources:
            n, e = source.outline(cs='latlon').T
            e, n = map(e,n)
            ax.fill(e, n, color=(0.5, 0.5, 0.5), lw=2)
            #add drawing of absolute nucleation point

    if syn_in.source() == 'DCSource':
            sources.append(DCSource(
                lat=float(syn_in.lat_0()),
                lon=float(syn_in.lon_0()),
                east_shift=float(syn_in.east_shift_0())*1000.,
                north_shift=float(syn_in.north_shift_0())*1000.,
                depth=syn_in.depth_syn_0()*1000.,
                strike=syn_in.strike_0(),
                dip=syn_in.dip_0(),
                rake=syn_in.rake_0(),
                time=util.str_to_time(syn_in.time_0()),
                magnitude=syn_in.magnitude_0()))
            for source in sources:
                ev = source.pyrocko_event()
                e, n = map(ev.lon,ev.lat)
                map.scatter(e, n, s=scale/200)


def load(filter, step=None, step_boot=None, booting_load=False):
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
            boot = False
            evpath = 'events/'+ str(sys.argv[1])
            C  = config.Config (evpath)
            Config = C.parseConfig('config')
            cfg = ConfigObj (dict=Config)
            n_bootstrap = cfg.UInt('n_bootstrap')
            dimx = int(Config['dimx'])
            dimy = int(Config['dimy'])
            data_int = None
            data = None
            data_boot = None
            data_int_boot = []
            datamax = 0
            phase = "P"
            for argv in sys.argv:
                if argv == "--phases:S":
                    phase = "S"
                if argv == "--phases:all":
                    phase = ""
                if argv == "--phases:P,S":
                    phase = ""
                if argv == "--phases:P":
                    phase = "P"
            if step is None:
                try:
                    pathlist = Path(rel).glob('%s-'+ str(sys.argv[5])+'*.ASC' % filter)
                except:
                    pathlist = Path(rel).glob('%s-*%s.ASC' % (filter,phase))
            else:

                try:
                    pathlist = Path(rel).glob('%s-'+ str(sys.argv[5])+'00%s_*.ASC' % (filter, step))
                except:
                    pathlist = Path(rel).glob('%s-*00%s_*%s.ASC' % (filter, step, phase))
            if booting_load is True:
                    pathlist = Path(rel).glob('%s-*00%s_*%s.ASC' % (filter, 0, phase))
            maxs = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    maxd = np.max(data[:, 2])
                    if maxs < maxd:
                        maxs = maxd
                        datamax = data[:, 2]
            if sys.argv[3] == 'max':
                if step is None:
                    try:
                        pathlist = Path(rel).glob('%s-'+ str(sys.argv[5])+'*.ASC' % filter)
                    except:
                        pathlist = Path(rel).glob('%s-*%s.ASC' % (filter, phase))
                else:
                    try:
                        pathlist = Path(rel).glob('%s-'+ str(sys.argv[5])+'00%s_*.ASC' % (filter, step))
                    except:
                        pathlist = Path(rel).glob('%s-*00%s_*%s.ASC' % (filter, step, phase))
                data_int = num.zeros(num.shape(data[:, 2]))
                for path in sorted(pathlist):
                        path_in_str = str(path)
                        data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                        i = 0
                        for k in np.nan_to_num(data[:,2]):
                            if k>data_int[i]:
                                data_int[i]= k
                            if num.max(datamax) == 0:
                                data_int[i]= 0
                            i = i+1
                try:
                    if booting_load is True:
                        boot = True
                        if step is None and step_boot is None:
                            pathlist = Path(rel).glob(('%s-*boot*'+'%s.ASC') % (filter, phase))
                        elif step is None:
                            pathlist = Path(rel).glob(('%s-*boot%s_*'+'%s.ASC') % (filter, step_boot, phase))
                        else:
                            pathlist = Path(rel).glob(('%s-*boot%s_*%s_'+'*%s.ASC') % (filter, step_boot, step, phase))
                        data_int_boot = num.zeros(num.shape(data[:, 2]))
                        for path in sorted(pathlist):
                                path_in_str = str(path)
                                data_boot = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                                i = 0
                                for k in np.nan_to_num(data_boot[:,2]):
                                    if k>data_int_boot[i]:
                                        data_int_boot[i]= k
                                #    if num.max(datamax) == 0:
                                #        data_int_boot[i]= 0
                                    i = i+1
                except IndexError:
                    pass

            if sys.argv[3] == 'combined':
                if step is None:
                    try:
                        pathlist = Path(rel).glob('%s-'+ str(sys.argv[5])+'*.ASC' % filter)
                    except:
                        pathlist = Path(rel).glob('%s-*%s.ASC' % (filter, phase))
                else:
                    try:
                        pathlist = Path(rel).glob('%s-'+ str(sys.argv[5])+'00%s_*.ASC' % (filter, step))
                    except:
                        pathlist = Path(rel).glob('%s-*00%s_*%s.ASC' % (filter, step, phase))
                data_int = num.zeros(num.shape(data[:, 2]))
                data_int_boot = num.ones(num.shape(data[:, 2]))

                for path in sorted(pathlist):
                        path_in_str = str(path)

                        if path_in_str[-14] is not "o":
                            data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                            i = 0
                            for k in np.nan_to_num(data[:,2]):
                                data_int_boot[i] = k+data_int_boot[i]
                                i = i+1
                try:
                    if booting_load is True:
                        if step is None and step_boot is None:
                            pathlist = Path(rel).glob(('%s-*boot*'+'%s.ASC') % (filter, phase))
                        elif step is None:
                            pathlist = Path(rel).glob(('%s-*boot%s_*'+'%s.ASC') % (filter, step_boot, phase))
                        else:
                            pathlist = Path(rel).glob(('%s-*boot%s_*%s_'+'%s.ASC') % (filter, step_boot, step, phase))
                        data_int_boot = num.ones(num.shape(data[:, 2]))
                        for path in sorted(pathlist):
                                path_in_str = str(path)
                                i = 0
                                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                                for k in np.nan_to_num(data[:,2]):
                                    data_int_boot[i] = k+data_int_boot[i]
                                    i = i+1
                    data_int = data_int_boot
                except IndexError:
                    pass

            return data, data_int, data_boot, data_int_boot, path_in_str, maxs, datamax

def make_map(data):
        eastings = data[:, 1]
        northings = data[:, 0]
        map = Basemap(projection='merc', llcrnrlon=num.min(eastings),
                      llcrnrlat=num.min(northings),
                      urcrnrlon=num.max(eastings),
                      urcrnrlat=num.max(northings),
                      resolution='h', epsg=3395)
        ratio_lat = num.max(northings)/num.min(northings)
        ratio_lon = num.max(eastings)/num.min(eastings)

        map.drawmapscale(num.min(eastings)+ratio_lon*0.25,
                         num.min(northings)+ratio_lat*0.25,
                         num.mean(eastings), num.mean(northings), 50)

        parallels = np.arange(num.min(northings),num.max(northings), int(ratio_lat))
        meridians = np.arange(num.min(eastings),num.max(eastings), int(ratio_lon))


        eastings, northings = map(eastings, northings)
        map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
        map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22, rotation=45)
        x, y = map(data[:,1], data[:,0])
        return map, x, y


def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.

    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")

    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)

    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):

        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu

    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi

    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi

    return (glon2, glat2, baz)

def great(m, startlon, startlat, azimuth,*args, **kwargs):
    glon1 = startlon
    glat1 = startlat
    glon2 = glon1
    glat2 = glat1

    step = 50

    glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    if azimuth-180 >= 0:
        while glon2 <= startlon:
            m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,**kwargs)
            azimuth = baz + 180.
            glat1, glon1 = (glat2, glon2)

            glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    elif azimuth-180 < 0:
        while glon2 >= startlon:
            m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,**kwargs)
            azimuth = baz + 180.
            glat1, glon1 = (glat2, glon2)

            glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)

def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])

    X,Y = m(X,Y)
    plt.plot(X,Y,color='gray',**kwargs)


def distance_time():

    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
    evpath = 'events/'+ str(sys.argv[1])

    C  = config.Config(evpath)
    Config = C.parseConfig('config')
    cfg = ConfigObj(dict=Config)
    step = cfg.UInt('step')
    step2 = cfg.UInt('step_f2')
    winlen = cfg.UInt('winlen')
    winlen2= cfg.UInt('winlen_f2')

    desired=[3,4]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
    desired=[7,8,9]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
    lat_ev, lon_ev = event_cor[1][0], event_cor[0][0]
    pathlist = Path(rel).glob('0-*.ASC')
    maxs = 0.
    if sys.argv[3] == 'combined':

        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                max = np.max(data[:, 2])
                if maxs < max:
                    maxs = max
                    datamax = np.max(data[:, 2])

        rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
        pathlist = Path(rel).glob('0-*.ASC')
        maxs = 0.
        data_int = num.zeros(num.shape(data[:, 2]))
        time_grid = num.zeros(num.shape(data[:, 2]))

        azis = []
        distances = []
        times = []
        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                data_int = data_int + data[:,2]
                for i in range(0, len(data[:, 2])):
                    if data_int[i] > datamax*0.1:
                        time_grid[i] = float(path_in_str[-8:-6]) * step
        for i in range(0, len(data[:, 2])):
            if data_int[i] > datamax*0.1:
                lats = data[i, 1]
                lons = data[i, 0]
                dist = orthodrome.distance_accurate50m(lats, lons,
                                                       lat_ev,
                                                       lon_ev)
                azis.append(toAzimuth(lat_ev, lon_ev,
                                      lats, lons))
                distances.append(dist)
                time = time_grid[i]
                times.append(time)
        datas = data_int
    if sys.argv[3] == 'max':

        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                max = np.max(data[:, 2])
                if maxs < max:
                    maxs = max
                    datamax = np.max(data[:, 2])

        rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
        pathlist = Path(rel).glob('0-*.ASC')
        maxs = 0.
        datas = []
        azis = []
        distances = []
        times = []
        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                for i in range(0, len(data[:, 2])):
                    if data[i, 2] > datamax*0.9:
                        lats = data[i, 1]
                        lons = data[i, 0]
                        datas.append(data[i, 2])
                        dist = orthodrome.distance_accurate50m(lats, lons,
                                                               lat_ev,
                                                               lon_ev)
                        azis.append(toAzimuth(lat_ev, lon_ev,
                                              lats, lons))
                        distances.append(dist)
                        time = float(path_in_str[-8:-6]) * step
                        times.append(time)

    if sys.argv[3] == 'stepwise':
        datamax = []
        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                max = np.max(data[:, 2])
                datamax.append(np.max(data[:, 2]))
        rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
        pathlist = Path(rel).glob('0-*.ASC')
        maxs = 0.
        datas = []
        azis = []
        distances = []
        times = []
        k = 0
        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                for i in range(0, len(data[:, 2])):
                    if data[i, 2] == datamax[k]:
                        lats = data[i, 1]
                        lons = data[i, 0]
                        datas.append(data[i, 2])
                        dist = orthodrome.distance_accurate50m(lats, lons,
                                                               lat_ev,
                                                               lon_ev)
                        azis.append(toAzimuth(lat_ev, lon_ev,
                                              lats, lons))
                        distances.append(dist)
                        time = float(path_in_str[-8:-6]) * step
                        times.append(time)
                k = k+1
        print(((distances[0]+distances[2])/2)/num.mean(time))

    if sys.argv[3] == 'stepwise_max':
        datamax = []
        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                max = np.max(data[:, 2])
                datamax.append(np.max(data[:, 2]))
        rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
        pathlist = Path(rel).glob('0-*.ASC')
        maxs = 0.
        datas = []
        azis = []
        distances = []
        times = []
        k = 0
        for path in sorted(pathlist):
                counter = 0
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                for i in range(0, len(data[:, 2])):
                    for kl in data[:, 2]:
                        if kl == datamax[k]:
                            counter = counter+1
                    if data[i, 2] == datamax[k]:
                        lats_list = []
                        lons_list = []
                        times_list = []
                        if counter == 0:
                            lats = data[i, 1]
                            lons = data[i, 0]
                            datas.append(data[i, 2])
                            dist = orthodrome.distance_accurate50m(lats, lons,
                                                                   lat_ev,
                                                                   lon_ev)
                            azis.append(toAzimuth(lat_ev, lon_ev,
                                                  lats, lons))
                            distances.append(dist)
                            time = float(path_in_str[-8:-6]) * step
                            times.append(time)
                        else:
                            lats_list.append(data[i, 1])
                            lons_list.append(data[i, 0])

                if counter != 0:
                    dist = orthodrome.distance_accurate50m(num.mean(lats_list), num.mean(lons_list),
                                                           lat_ev,
                                                           lon_ev)
                    distances.append(dist)

                    time = float(path_in_str[-8:-6]) * step
                    times.append(time)
                k = k+1
        print(num.mean(distances)/num.mean(time))

    plt.figure()
    plt.scatter(distances, times, s=100)
    plt.show()

    plt.figure()
    plt.scatter(azis, times, s=100)
    plt.show()


def distance_time_bootstrap():

    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
    evpath = 'events/'+ str(sys.argv[1])
    C  = config.Config(evpath)
    Config = C.parseConfig('config')
    cfg = ConfigObj(dict=Config)
    step = cfg.UInt('step')
    step2 = cfg.UInt('step_f2')
    winlen = cfg.UInt('winlen')
    winlen2= cfg.UInt('winlen_f2')
    n_bootstrap = cfg.UInt('n_bootstrap')

    desired=[3,4]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
    desired=[7,8,9]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
    lat_ev, lon_ev = event_cor[1][0], event_cor[0][0]
    pathlist_main = Path(rel).glob('0-*.ASC')

    maxs = 0.
    if sys.argv[3] == 'combined':

        for path in sorted(pathlist_main):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                max = np.max(data[:, 2])
                if maxs < max:
                    maxs = max
                    datamax = np.max(data[:, 2])

        rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
        pathlist_main = Path(rel).glob('0-*.ASC')
        data_int_boot = []
        maxs = 0.
        data_int = num.zeros(num.shape(data[:, 2]))
        time_grid_boot = []

        for n in range(0, n_bootstrap+1):
            data_int_boot.append(num.zeros(num.shape(data[:, 2])))
            time_grid_boot.append(num.zeros(num.shape(data[:, 2])))

        time_grid = num.zeros(num.shape(data[:, 2]))
        for path in sorted(pathlist_main):
                path_in_str = str(path)
                if path_in_str[-14] is not "o":
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

        azis = []
        distances = []
        times = []
        list_boots = num.arange(0, n_bootstrap)

        pathlist_main = Path(rel).glob('0-*.ASC')

        for path in sorted(pathlist_main):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)

                if path_in_str[-14] is not "o":
                    data_int = data_int + data[:,2]
                    for i in range(0, len(data[:, 2])):
                        if data_int[i] > datamax*0.1:
                            time_grid[i] = float(path_in_str[-8:-6]) * step

                else:
                    bnr = int(path_in_str[-11])
                    data_int_boot[bnr] = data_int_boot[bnr] + data[:,2]
                    for i in range(0, len(data[:, 2])):
                        if data_int_boot[bnr][i] > datamax*0.1:
                            time_grid_boot[bnr][i] = float(path_in_str[-8:-6]) * step

        for i in range(0, len(data[:, 2])):
            if data_int[i] > datamax*0.1:
                lats = data[i, 1]
                lons = data[i, 0]
                dist = orthodrome.distance_accurate50m(lats, lons,
                                                       lat_ev,
                                                       lon_ev)
                azis.append(toAzimuth(lat_ev, lon_ev,
                                      lats, lons))
                distances.append(dist)
                time = time_grid[i]
                times.append(time)

        datas = data_int
        distances_boot = []
        times_boot = []
        for bnr in range(0, n_bootstrap+1):
            time_boot = []
            dist_boot = []
            for i in range(0, len(data[:, 2])):
                if data_int_boot[bnr][i] > datamax*0.1:
                    lats = data[i, 1]
                    lons = data[i, 0]
                    dist = orthodrome.distance_accurate50m(lats, lons,
                                                           lat_ev,
                                                           lon_ev)
                    azis.append(toAzimuth(lat_ev, lon_ev,
                                          lats, lons))
                    dist_boot.append(dist)
                    time = time_grid[i]
                    time_boot.append(time)
            times_boot.append(time_boot)
            distances_boot.append(dist_boot)
    if sys.argv[3] == 'max':

        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                max = np.max(data[:, 2])
                if maxs < max:
                    maxs = max
                    datamax = np.max(data[:, 2])

        rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
        pathlist = Path(rel).glob('0-*.ASC')
        maxs = 0.
        datas = []
        azis = []
        distances = []
        times = []
        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                for i in range(0, len(data[:, 2])):
                    if data[i, 2] > datamax*0.9:
                        lats = data[i, 1]
                        lons = data[i, 0]
                        datas.append(data[i, 2])
                        dist = orthodrome.distance_accurate50m(lats, lons,
                                                               lat_ev,
                                                               lon_ev)
                        azis.append(toAzimuth(lat_ev, lon_ev,
                                              lats, lons))
                        distances.append(dist)
                        time = float(path_in_str[-8:-6]) * step
                        times.append(time)

    if sys.argv[3] == 'stepwise':
        datamax = []
        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                max = np.max(data[:, 2])
                datamax.append(np.max(data[:, 2]))
        rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
        pathlist = Path(rel).glob('0-*.ASC')
        maxs = 0.
        datas = []
        azis = []
        distances = []
        times = []
        k = 0
        for path in sorted(pathlist):
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                for i in range(0, len(data[:, 2])):
                    if data[i, 2] == datamax[k]:
                        lats = data[i, 1]
                        lons = data[i, 0]
                        datas.append(data[i, 2])
                        dist = orthodrome.distance_accurate50m(lats, lons,
                                                               lat_ev,
                                                               lon_ev)
                        azis.append(toAzimuth(lat_ev, lon_ev,
                                              lats, lons))
                        distances.append(dist)
                        time = float(path_in_str[-8:-6]) * step
                        times.append(time)

                k = k+1
    plt.figure()
    colors = iter(cm.rainbow(np.linspace(0, 1, n_bootstrap+1)))

    plt.scatter(distances, times, c='k', s=100)
    for n in range(0, n_bootstrap+1):
        plt.scatter(distances_boot[n], times_boot[n], c=next(colors), s=100)

    plt.show()

def plot_cluster():

    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
    desired=[3,4]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
    desired=[7,8,9]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
    map = Basemap(width=21000000,height=21000000,
                resolution='l',projection='aeqd',\
                lat_ts=event_cor[0][0],lat_0=event_cor[0][0],lon_0=event_cor[1][0])
    map.fillcontinents(zorder=-1)
    map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=[0,0,0,1])
    x, y = map(event_cor[1][0],event_cor[0][0])
    ax = plt.gca()
    np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
    beach1 = beach(np1, xy=(x, y), width=900030)
    ax.add_collection(beach1)
    pathlist = Path(rel).glob('*.dat')
    i=0
    for path in sorted(pathlist):
        path_in_str = str(path)
        i = i+1
    colors = iter(cm.rainbow(np.linspace(0, 1, i)))
    pathlist = Path(rel).glob('*.dat')
    for path in sorted(pathlist):
        path_in_str = str(path)
        data = num.loadtxt(path_in_str, delimiter=' ', usecols=(0,3,4))
        try:
            lons = data[:,2]
            lats = data[:,1]

        except:
            lons = data[2]
            lats = data[1]

        x, y = map(lons,lats)
        map.scatter(x,y,30,marker='o',c=next(colors))
        try:
            plt.text(x[0],y[0],'r'+str(data[0,0])[:], fontsize=12)
        except:
            plt.text(x,y,'r'+str(data[0])[0:2], fontsize=12)
            pass
        lon_0, lat_0 = event_cor[1][0],event_cor[0][0]
        x,y=map(lon_0,lat_0)
        degree_sign= u'\N{DEGREE SIGN}'
        x2,y2 = map(lon_0,lat_0-20)
        plt.text(x2,y2,'20'+degree_sign, fontsize=22,color='blue')
        circle1 = plt.Circle((x, y), y2-y, color='blue',fill=False, linestyle='dashed')
        ax.add_patch(circle1)
        x,y=map(lon_0,lat_0)
        x2,y2 = map(lon_0,lat_0-60)
        plt.text(x2,y2,'60'+degree_sign, fontsize=22,color='blue')
        circle2 = plt.Circle((x, y), y2-y, color='blue',fill=False, linestyle='dashed')
        ax.add_patch(circle2)
        x,y=map(lon_0,lat_0)
        x2,y2 = map(lon_0,lat_0-90)
        plt.text(x2,y2,'90'+degree_sign, fontsize=22,color='blue')
        circle2 = plt.Circle((x, y), y2-y, color='blue',fill=False, linestyle='dashed')
        ax.add_patch(circle2)
        x,y=map(lon_0,lat_0)
        x2,y2 = map(lon_0,lat_0-94)
        circle2 = plt.Circle((x, y), y2-y, color='red',fill=False, linestyle='dashed')
        ax.add_patch(circle2)
        x,y=map(lon_0,lat_0)
        x2,y2 = map(lon_0,lat_0-22)
        circle2 = plt.Circle((x, y), y2-y, color='red',fill=False, linestyle='dashed')
        ax.add_patch(circle2)
    plt.show()


def plot_timeshift_map():

    import _pickle as pickle

    evpath = 'events/'+ str(sys.argv[1])
    C = config.Config(evpath)
    Config = C.parseConfig('config')
    cfg = ConfigObj(dict=Config)
    sembpath = evpath + '/work/semblance'
    stations = []
    refs = []
    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

    if cfg.Bool('synthetic_test') is True:
        Syn_in = C.parseConfig('syn')
        syn_in = SynthCfg(Syn_in)
        lat_ev = float(syn_in.lat_0())
        lon_ev = float(syn_in.lon_0())
    else:
        event = 'events/' + str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
        desired=[3,4]
        with open(event, 'r') as fin:
            reader = csv.reader(fin)
            event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
            lat_ev, lon_ev = event_cor[1][0], event_cor[0][0]


    pathlist = Path(rel).glob('*.shift*l0*')
    for path in sorted(pathlist):
            path_in_str = str(path)
            if path_in_str[-1] != "s":
                f = open(path_in_str, 'rb')
                refshifts = pickle.load(f)
                f.close()
                for s in refshifts.values():
                    refs.append(s)

            else:
                f = open(path_in_str, 'rb')
                refshifts_stations = pickle.load(f)
                f.close()
                for s in refshifts_stations.values():
                    stations.append(s)

    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
    desired=[3,4]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
    desired=[7,8,9]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
    map = Basemap(width=21000000,height=21000000,
                resolution='l',projection='aeqd',\
                lat_ts=event_cor[0][0],lat_0=event_cor[0][0],lon_0=event_cor[1][0])
    map.fillcontinents(zorder=-1)
    map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=[0,0,0,1])
    x, y = map(event_cor[1][0],event_cor[0][0])
    ax = plt.gca()
    np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
    beach1 = beach(np1, xy=(x, y), width=900030)
    ax.add_collection(beach1)
    pathlist = Path(rel).glob('*.dat')
    i=0
    minima = min(refs)
    maxima = max(refs)
    import matplotlib

    norm = matplotlib.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.jet)
    for st, ref in zip(stations, refs):


        x, y = map(st[1],st[0])

        map.scatter(x,y,30,marker='o',c=mapper.to_rgba(ref))

    lon_0, lat_0 = event_cor[1][0],event_cor[0][0]
    x,y=map(lon_0,lat_0)
    degree_sign= u'\N{DEGREE SIGN}'
    x2,y2 = map(lon_0,lat_0-20)
    plt.text(x2,y2,'20'+degree_sign, fontsize=22,color='blue')
    circle1 = plt.Circle((x, y), y2-y, color='blue',fill=False, linestyle='dashed')
    ax.add_patch(circle1)
    x,y=map(lon_0,lat_0)
    x2,y2 = map(lon_0,lat_0-60)
    plt.text(x2,y2,'60'+degree_sign, fontsize=22,color='blue')
    circle2 = plt.Circle((x, y), y2-y, color='blue',fill=False, linestyle='dashed')
    ax.add_patch(circle2)
    x,y=map(lon_0,lat_0)
    x2,y2 = map(lon_0,lat_0-90)
    plt.text(x2,y2,'90'+degree_sign, fontsize=22,color='blue')
    circle2 = plt.Circle((x, y), y2-y, color='blue',fill=False, linestyle='dashed')
    ax.add_patch(circle2)
    x,y=map(lon_0,lat_0)
    x2,y2 = map(lon_0,lat_0-94)
    circle2 = plt.Circle((x, y), y2-y, color='red',fill=False, linestyle='dashed')
    ax.add_patch(circle2)
    x,y=map(lon_0,lat_0)
    x2,y2 = map(lon_0,lat_0-22)
    circle2 = plt.Circle((x, y), y2-y, color='red',fill=False, linestyle='dashed')
    ax.add_patch(circle2)
    plt.show()




    pathlist = Path(rel).glob('*.shift*h1*')
    for path in sorted(pathlist):
            path_in_str = str(path)
            if path_in_str[-1] != "s":
                f = open(path_in_str, 'rb')
                refshifts = pickle.load(f)
                f.close()
                for s in refshifts.values():
                    refs.append(s)

            else:
                f = open(path_in_str, 'rb')
                refshifts_stations = pickle.load(f)
                f.close()
                for s in refshifts_stations.values():
                    stations.append(s)

    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
    desired=[3,4]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
    desired=[7,8,9]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
    map = Basemap(width=21000000,height=21000000,
                resolution='l',projection='aeqd',\
                lat_ts=event_cor[0][0],lat_0=event_cor[0][0],lon_0=event_cor[1][0])
    map.fillcontinents(zorder=-1)
    map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=[0,0,0,1])
    x, y = map(event_cor[1][0],event_cor[0][0])
    ax = plt.gca()
    np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
    beach1 = beach(np1, xy=(x, y), width=900030)
    ax.add_collection(beach1)
    pathlist = Path(rel).glob('*.dat')
    i=0
    minima = min(refs)
    maxima = max(refs)
    import matplotlib

    norm = matplotlib.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.jet)
    for st, ref in zip(stations, refs):


        x, y = map(st[1],st[0])

        map.scatter(x,y,30,marker='o',c=mapper.to_rgba(ref))

    lon_0, lat_0 = event_cor[1][0],event_cor[0][0]
    x,y=map(lon_0,lat_0)
    degree_sign= u'\N{DEGREE SIGN}'
    x2,y2 = map(lon_0,lat_0-20)
    plt.text(x2,y2,'20'+degree_sign, fontsize=22,color='blue')
    circle1 = plt.Circle((x, y), y2-y, color='blue',fill=False, linestyle='dashed')
    ax.add_patch(circle1)
    x,y=map(lon_0,lat_0)
    x2,y2 = map(lon_0,lat_0-60)
    plt.text(x2,y2,'60'+degree_sign, fontsize=22,color='blue')
    circle2 = plt.Circle((x, y), y2-y, color='blue',fill=False, linestyle='dashed')
    ax.add_patch(circle2)
    x,y=map(lon_0,lat_0)
    x2,y2 = map(lon_0,lat_0-90)
    plt.text(x2,y2,'90'+degree_sign, fontsize=22,color='blue')
    circle2 = plt.Circle((x, y), y2-y, color='blue',fill=False, linestyle='dashed')
    ax.add_patch(circle2)
    x,y=map(lon_0,lat_0)
    x2,y2 = map(lon_0,lat_0-94)
    circle2 = plt.Circle((x, y), y2-y, color='red',fill=False, linestyle='dashed')
    ax.add_patch(circle2)
    x,y=map(lon_0,lat_0)
    x2,y2 = map(lon_0,lat_0-22)
    circle2 = plt.Circle((x, y), y2-y, color='red',fill=False, linestyle='dashed')
    ax.add_patch(circle2)
    plt.show()

def plot_movie():

    evpath = 'events/'+ str(sys.argv[1])
    C  = config.Config (evpath)
    Config = C.parseConfig ('config')
    cfg = ConfigObj (dict=Config)
    step = cfg.UInt ('step')
    step2 = cfg.UInt ('step_f2')
    winlen = cfg.UInt ('winlen')
    winlen2= cfg.UInt ('winlen_f2')
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
            pathlist = Path(rel).glob('0-*.ASC')
            maxs = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]
            pathlist = Path(rel).glob('0-*.ASC')
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    eastings = data[:,1]
                    northings =  data[:,0]
                    plt.figure()
                    time = float(path_in_str[-8:-6])* step

                    map = Basemap(projection='merc', llcrnrlon=num.min(eastings),
                                  llcrnrlat=num.min(northings),
                                  urcrnrlon=num.max(eastings),
                                  urcrnrlat=num.max(northings),
                                  resolution='h', epsg=3395)
                    ratio_lat = num.max(northings)/num.min(northings)
                    ratio_lon = num.max(eastings)/num.min(eastings)

                    map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)

                    parallels = np.arange(num.min(northings),num.max(northings),0.2)
                    meridians = np.arange(num.min(eastings),num.max(eastings),0.2)
                    xpixels = 1000
                    map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
                    eastings, northings = map(eastings, northings)
                    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
                    map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
                    x, y = map(data[::10,1], data[::10,0])
                    mins = np.max(data[:,2])
                    plt.tricontourf(x,y, data[::10,2], cmap='hot', vmin=0., vmax=maxs)
                    plt.colorbar()
                    plt.title(path_in_str+'first filter')
                    plt.savefig('time:'+str(time)+'_f1'+'.png', bbox_inches='tight')
                    plt.close()
            try:
                pathlist = Path(rel).glob('1-*.ASC')
                for path in sorted(pathlist):
                        path_in_str = str(path)
                        time = float(path_in_str[-8:-6])* step
                        data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                        eastings = data[:,1]
                        northings =  data[:,0]
                        plt.figure()

                        map = Basemap(projection='merc', llcrnrlon=num.min(eastings),
                                      llcrnrlat=num.min(northings),
                                      urcrnrlon=num.max(eastings),
                                      urcrnrlat=num.max(northings),
                                      resolution='h', epsg=3395)
                        ratio_lat = num.max(northings)/num.min(northings)
                        ratio_lon = num.max(eastings)/num.min(eastings)

                        map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)
                        parallels = np.arange(num.min(northings),num.max(northings),0.2)
                        meridians = np.arange(num.min(eastings),num.max(eastings),0.2)
                        xpixels = 1000
                        map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
                        eastings, northings = map(eastings, northings)
                        map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
                        map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
                        x, y = map(data[::10,1], data[::10,0])
                        mins = np.max(data[:,2])
                        plt.tricontourf(x,y, data[::10,2], cmap='hot', vmin=0., vmax=maxs)
                        plt.colorbar()
                        plt.title(path_in_str+'second filter')
                        plt.savefig('time:'+str(time)+'_f2'+'.png', bbox_inches='tight')
                        plt.close()
            except:
                pass

        else:
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/' + str(sys.argv[3])
            pathlist = Path(rel).glob('**/*.ASC')
            for path in sorted(pathlist):
                try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    eastings = data[:,1]
                    northings =  data[:,0]
                    plt.figure()
                    map = Basemap(projection='merc', llcrnrlon=num.min(eastings),
                                  llcrnrlat=num.min(northings),
                                  urcrnrlon=num.max(eastings),
                                  urcrnrlat=num.max(northings),
                                  resolution='h', epsg=3395)
                    ratio_lat = num.max(northings)/num.min(northings)
                    ratio_lon = num.max(eastings)/num.min(eastings)

                    map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)
                    parallels = np.arange(num.min(northings),num.max(northings),0.2)
                    meridians = np.arange(num.min(eastings),num.max(eastings),0.2)

                    eastings, northings = map(eastings, northings)
                    map.drawcoastlines(color='b',linewidth=3)
                    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
                    map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
                    x, y = map(data[::5,1], data[::5,0])
                    mins = np.max(data[:,3])
                    plt.tricontourf(x,y, data[::5,3], vmin=mins*0.6)
                    plt.title(path_in_str)
                    plt.savefig(path_in_str+'_f1'+'.pdf', bbox_inches='tight')
                    plt.close()
                except:
                    plt.close()
                    pass


def integrated_scatter():


    evpath = 'events/'+ str(sys.argv[1])
    C  = config.Config (evpath)
    Config = C.parseConfig ('config')
    cfg = ConfigObj (dict=Config)
    step = cfg.Float ('step')
    step2 = cfg.Float ('step_f2')

    if len(sys.argv)<5:
        print("missing input arrayname and or depth")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
            import matplotlib
            matplotlib.rcParams.update({'font.size': 32})
            pathlist = Path(rel).glob('0-'+'*boot*.ASC')
            maxs = 0.
            counter = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    counter =+ 1
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]

            pathlist = Path(rel).glob('0-'+'*boot*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            data_old = num.zeros(num.shape(data[:, 2]))
            time_grid = num.zeros(num.shape(data[:, 2]))
            times = []

            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])
                    time = float(path_in_str[-8:-6])* step
                    times.append(time)
                    for i in range(0, num.shape(data[:, 2])[0]):
                        if data[i,2] > data_old[i] and time_grid[i] == 0:
                            time_grid[i] = time
                            data_old[i] = data[i,2]

            eastings = data[:,1]
            northings = data[:,0]
            xpixels = 1000

            plt.figure()

            map, x, y = make_map(data)
            mins = np.max(data[:,2])
            size =(data_int/np.max(data_int))*300
            times_idx = np.where(time_grid==0)
            time_grid[times_idx] = 'NaN'
            ps = map.scatter(x,y,marker='o',c=time_grid, s=size, cmap='jet')
            cb = plt.colorbar(orientation="horizontal")
            cb.outline.set_visible(False)
            cb.set_label('Time ->',fontsize=22)
            plt.title(path_in_str)

            xpixels = 6000
            eastings = data[:,1]
            northings = data[:,0]
            map.arcgisimage(service='World_Shaded_Relief',
                            xpixels=xpixels, verbose=False)
            plt.show()
            pathlist = Path(rel).glob('1-'+'*boot*.ASC')
            maxs = 0.
            counter = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    counter =+ 1
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]

            pathlist = Path(rel).glob('1-'+'*boot*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            data_old = num.zeros(num.shape(data[:, 2]))
            time_grid = num.zeros(num.shape(data[:, 2]))
            times = []

            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])
                    time = float(path_in_str[-8:-6])* step
                    times.append(time)
                    for i in range(0, num.shape(data[:, 2])[0]):
                        if data[i,2] > data_old[i]:
                            time_grid[i] = time
                            data_old[i] = data[i,2]

            eastings = data[:,1]
            northings = data[:,0]
            xpixels = 1000

            plt.figure()

            map, x, y = make_map(data)
            mins = np.max(data[:,2])
            size =(data_int/np.max(data_int))*300
            times_idx = np.where(time_grid==0)
            time_grid[times_idx] = 'NaN'
            ps = map.scatter(x,y,marker='o',c=time_grid, s=size, cmap='jet')
            cb = plt.colorbar(orientation="horizontal")
            cb.outline.set_visible(False)
            cb.set_label('Time ->',fontsize=22)
            plt.title(path_in_str)

            xpixels = 6000
            eastings = data[:,1]
            northings = data[:,0]
            map.arcgisimage(service='World_Shaded_Relief',
                            xpixels=xpixels, verbose=False)
            plt.show()

def beampower():
        rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
        pathlist = Path(rel).glob('r*/beam.mseed')
        for path in sorted(pathlist):
                path_in_str = str(path)
                tr_bp = io.load(path_in_str)[0]
                tr_bp.ydata = tr_bp.ydata*0.
        pathlist = Path(rel).glob('r*/beam.mseed')
        for path in sorted(pathlist):
                path_in_str = str(path)
                tr = io.load(path_in_str)[0]
                tr.ydata = abs(tr.ydata)
                tr_bp.add(tr)
        trace.snuffle(tr_bp)
        bp_diff_tr = tr_bp.copy()
        bp_diff_tr.ydata = num.diff(tr_bp.ydata)
        trace.snuffle(tr_bp_diff)

def spec(tr):
        f, a = tr.spectrum(pad_to_pow2=True)
        return (f, a)

def inspect_spectrum():
        from pyrocko import cake
        event = model.load_events('events/'+ str(sys.argv[1]) + '/data/event.pf')[0]
        rel = 'events/'+ str(sys.argv[1]) + '/data/'
        traces = io.load(rel+'traces_velocity.mseed')
        stations = model.load_stations(rel+'stations.txt')
        earth = cake.load_model('ak135-f-continental.m')
        for tr in traces:
            for st in stations:
                if tr.station == st.station and tr.location == st.location:
                        distances = [st.distance_to(event)* cake.m2d, st.distance_to(event)* cake.m2d]
                        Phase = cake.PhaseDef('P')
                        rays = earth.arrivals(
                            phases=Phase,
                            distances=distances,
                            zstart=event.depth*2,
                            zstop=0.0)
                        time = rays[0].t+event.time
                        tr.ydata = tr.ydata.astype(num.float)
                        tr.ydata -= tr.ydata.mean()
                        tr_spec, a = spec(tr)
                        tr.ydata = abs(a)

        trace.snuffle(traces)



def plot_integrated_movie():
    evpath = 'events/' + str(sys.argv[1])
    C = config.Config(evpath)
    Config = C.parseConfig('config')
    cfg = ConfigObj(dict=Config)
    step = cfg.UInt('step')
    step2 = cfg.UInt('step_f2')
    duration = cfg.UInt('duration')
    forerun = cfg.UInt('forerun')
    plt_time = False
    boot = False
    ntimes = int((forerun+duration)/step)

    dimx = int(Config['dimx'])
    dimy = int(Config['dimy'])
    filters = cfg.String('filters')
    filters = int(filters)
    max_rel = False
    try:
        for argv in sys.argv:
            if argv == 'boot':
                boot = True
            if argv == '--time':
                plt_time = True
            if argv == "--max_rel":
                max_rel = True
    except:
        pass
    if plt_time is True:
        fig = plt.figure()

    for filterindex in range(0, filters):
        data_all, data_int_all, data_boot, data_int_boot, path_in_str, maxs, datamax = load(filterindex)
        cmaps = ['Blues', 'Greens', 'Reds', 'Purples', 'Greys', 'Wistia', 'bone', 'copper', 'dusk']

        levels = np.linspace(num.min(data_int_all), maxs, 20)
        for i in range(0, ntimes):
            if len(sys.argv) < 4:
                print("missing input arrayname")
            else:
                    data, data_int, data_boot, data_int_boot, path_in_str, maxsb, datamaxb = load(filterindex, step=i)
                    if plt_time is False:
                        fig = plt.figure()
                    try:
                        ax = fig.axes[0]
                    except:
                        ax = fig.axes
                    map, x, y = make_map(data)
                    xmax = num.max(x)
                    xmin = num.min(x[num.nonzero(x)])
                    ymax = num.max(y)
                    ymin = num.min(y[num.nonzero(y)])
                    scale = (xmax/xmin)*(ymax/ymin)*10
                    triang = tri.Triangulation(x, y)
                    isbad = np.less(data_int, num.min(data_int))
                    mask = np.all(np.where(isbad[triang.triangles],
                                           True, False), axis=1)
                    triang.set_mask(mask)
                    if plt_time is False:
                        plt.tricontourf(triang, data_int, vmax=maxs, vmin=0,
                                        cmap=cm.coolwarm)
                        m = plt.cm.ScalarMappable(cmap=cm.coolwarm)
                        m.set_array(data_int)
                        m.set_clim(0., maxs)
                        plt.colorbar(m, orientation="horizontal",
                                     boundaries=np.linspace(0, maxs, 10))
                        plt.title(path_in_str)

                    else:
                        colors = []

                        plt.tricontourf(triang, data_int,
                                        cmap=cmaps[i], levels=levels)

                    event = 'events/' + str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
                    draw_beach(ax, scale, map, event)

                    if boot is True:
                        n_bootstrap = cfg.UInt('n_bootstrap')
                        for iboot in range(0, n_bootstrap):
                            datab, data_intb, data_boot, data_int_boot, path_in_strb, maxsb, datamaxb = load(filterindex, step=i, step_boot=iboot, booting_load=True)
                            where_are_NaNs = num.isnan(data_int_boot)
                            data_int_boot[where_are_NaNs] = 0
                            data_int_boot = num.reshape(data_int_boot, (dimx,
                                                                        dimy))
                            xc = num.reshape(x, (dimx, dimy))
                            yc = num.reshape(y, (dimx, dimy))
                            plot_comb_bs = False
                            plot_ind_bs = True
                            if plot_ind_bs is True:
                                try:
                                    if plt_time is False:
                                        cp= plt.contour(xc, yc, data_int_boot, levels=[num.std(data_int_boot)*2])
                                    else:
                                        cmap = cm.get_cmap(cmaps[i])
                                        cp= plt.contour(xc, yc, data_int_boot, levels=[num.std(data_int_boot)*2], colors=[cmap(1.0)])

                                except ValueError:
                                    pass

                        if plot_comb_bs is True:
                            datab, data_intb, data_boot, data_int_boot, path_in_strb, maxsb, datamaxb = load(filterindex)

                            offset = [num.mean(x), num.mean(y)]
                            data_int_boot = num.reshape(data_int_boot, (dimx,
                                                                        dimy))
                            xc = num.reshape(x, (dimx, dimy))
                            yc = num.reshape(y, (dimx, dimy))
                            cp= plt.contour(xc, yc, data_int_boot, levels=[num.std(data_intb)*2])
                            ax.clabel(cp, fmt='%2.1f', colors='w', fontsize=14)

                    for argv in sys.argv:
                        if argv == "--topography":
                            try:
                                xpixels = 1000
                                map.arcgisimage(service='World_Shaded_Relief',
                                                xpixels = xpixels,
                                                verbose= False)
                            except:
                                pass

                    if plt_time is False:
                        if cfg.Bool('synthetic_test') is True:
                            Syn_in = C.parseConfig('syn')
                            syn_in = SynthCfg(Syn_in)
                            ax = plt.gca()

                            draw_sources(ax, syn_in, map, scale)
                        plt.savefig('time:'+str(i)+'_f1'+'.png', bbox_inches='tight')
                        plt.show()

        if plt_time is True:
            if cfg.Bool('synthetic_test') is True:
                Syn_in = C.parseConfig('syn')
                syn_in = SynthCfg(Syn_in)
                draw_sources(ax, syn_in, map, scale)
            plt.show()


def plot_integrated():
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
            evpath = 'events/' + str(sys.argv[1])
            C = config.Config(evpath)
            Config = C.parseConfig('config')
            cfg = ConfigObj(dict=Config)
            filters = cfg.String('filters')
            filters = int(filters)
            boot = False
            try:
                for argv in sys.argv:
                    if argv == 'boot':
                        boot = True
            except:
                pass

            for filterindex in range(0, filters):
                data, data_int, data_boot, data_int_boot, path_in_str, maxs, datamax = load(filterindex)
                rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

                dimx = int(Config['dimx'])
                dimy = int(Config['dimy'])

                plt.figure()
                ax = plt.gca()

                map, x, y = make_map(data)
                xmax = num.max(x)
                xmin = num.min(x[num.nonzero(x)])
                ymax = num.max(y)
                ymin = num.min(y[num.nonzero(y)])
                scale = (xmax/xmin)*(ymax/ymin)*10

                triang = tri.Triangulation(x, y)
                isbad = np.less(data_int, num.max(data_int)*0.001)
                mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
                triang.set_mask(mask)
                plt.tricontourf(triang, data_int, cmap='cool')
                plt.colorbar(orientation="horizontal")
                plt.title(path_in_str)

                if boot is True:
                    plot_comb_bs = False
                    plot_ind_bs = False
                    print('plotting boot')
                    for argv in sys.argv:
                        if argv == "--induvidual":
                            plot_ind_bs = True
                        else:
                            plot_comb_bs = True
                    n_bootstrap = cfg.UInt('n_bootstrap')
                    for iboot in range(0, n_bootstrap):
                        datab, data_intb, data_boot, data_int_boot, path_in_strb, maxsb, datamaxb = load(filterindex, step_boot=iboot, booting_load=True)

                        where_are_NaNs = num.isnan(data_int_boot)
                        data_int_boot[where_are_NaNs] = 0
                        data_int_boot = num.reshape(data_int_boot, (dimx,
                                                                    dimy))
                        xc = num.reshape(x, (dimx, dimy))
                        yc = num.reshape(y, (dimx, dimy))

                        if plot_ind_bs is True:

                            try:
                                cp = plt.contour(xc, yc, data_int_boot, levels=[num.std(data_intb)*4])
                            except ValueError:
                                pass

                    if plot_comb_bs is True:
                        datab, data_intb, data_boot, data_int_boot, path_in_strb, maxsb, datamaxb = load(filterindex, booting_load=True)

                        offset = [num.mean(x), num.mean(y)]
                        data_int_boot = num.reshape(data_int_boot, (dimx,
                                                                    dimy))
                        xc = num.reshape(x, (dimx, dimy))
                        yc = num.reshape(y, (dimx, dimy))
                        cp = plt.contour(xc, yc, data_int_boot, levels=[num.std(data_intb)*2])
                        ax.clabel(cp, fmt='%2.1f', colors='w', fontsize=14)

                event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
                draw_beach(ax, scale, map, event)

                for argv in sys.argv:
                    if argv == "--topography":
                        try:
                            xpixels = 1000
                            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
                        except:
                            pass

                if cfg.Bool('synthetic_test') is True:
                    Syn_in = C.parseConfig('syn')
                    syn_in = SynthCfg(Syn_in)
                    draw_sources(ax, syn_in, map, scale)
                plt.show()


def plot_time():
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

            try:
                pathlist = Path(rel).glob('0-'+ str(sys.argv[5])+'.ASC')
            except:
                pathlist = Path(rel).glob('0-*.ASC')
            maxs = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]

            try:
                pathlist = Path(rel).glob('0-'+ str(sys.argv[5])+'.ASC')
            except:
                pathlist = Path(rel).glob('0-*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    i = 0
                    for k in np.nan_to_num(data[:,2]):
                        if k>data_int[i]:
                            data_int[i]= k
                        i = i+1

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()
            map, x, y = make_map(data)

            mins = np.max(data[:,2])
            triang = tri.Triangulation(x, y)
            isbad = np.less(data_int, 0.085)
            mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            levels = np.arange(0., 1.05, 0.025)
            triang.set_mask(mask)
            plt.tricontourf(triang, data_int, cmap='cool')
            plt.colorbar(orientation="horizontal")
            plt.title(path_in_str)
            event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
            desired=[3,4]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
            desired=[7,8,9]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
            x, y = map(event_cor[1][0],event_cor[0][0])
            ax = plt.gca()
            np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
            beach1 = beach(np1, xy=(x, y), width=0.09)
            ax.add_collection(beach1)
            for argv in sys.argv:
                if argv == "--topography":
                    try:
                        xpixels = 1000
                        map.arcgisimage(service='World_Shaded_Relief',
                                        xpixels = xpixels,
                                        verbose= False)
                    except:
                        pass

            plt.show()

            try:
                pathlist = Path(rel).glob('1-'+ str(sys.argv[5])+'.ASC')
            except:
                pathlist = Path(rel).glob('1-*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()
            map, x, y = make_map(data)

            mins = np.max(data[:,2])

            triang = tri.Triangulation(x, y)
            isbad = np.less(data_int, 0.01)
            mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang.set_mask(mask)
            plt.tricontourf(triang, data_int, cmap='YlOrRd')
            event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
            desired=[3,4]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
            desired=[7,8,9]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
            x, y = map(event_cor[1][0],event_cor[0][0])
            ax = plt.gca()
            np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
            beach1 = beach(np1, xy=(x, y), width=0.09)
            ax.add_collection(beach1)
            plt.colorbar()
            plt.title(path_in_str)
            for argv in sys.argv:
                if argv == "--topography":
                    try:
                        xpixels = 1000
                        map.arcgisimage(service='World_Shaded_Relief',
                                        xpixels = xpixels,
                                        verbose= False)
                    except:
                        pass

            plt.show()


def plot_semb_equal():
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

            try:
                pathlist = Path(rel).glob('0-'+ str(sys.argv[5])+'.ASC')
            except:
                pathlist = Path(rel).glob('0-*.ASC')
            maxs = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]

            try:
                pathlist = Path(rel).glob('0-'+ str(sys.argv[5])+'.ASC')
            except:
                pathlist = Path(rel).glob('0-*.ASC')


            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()
            map, x, y = make_map(data)

            mins = np.max(data[:,2])
            triang = tri.Triangulation(x, y)
            #isbad = np.less(data_int, 0.085)
            #mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            levels = np.arange(0., 1.05, 0.025)
            #triang.set_mask(mask)
            for path in sorted(pathlist):
                    data_int = num.zeros(num.shape(data[:, 2]))
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    i = 0
                    for k in np.nan_to_num(data[:,2]):
                        if k>data_int[i]:
                            data_int[i]= k
                        i = i+1
                    plt.tricontourf(triang, data_int, cmap='cool')
            plt.colorbar(orientation="horizontal")
            plt.title(path_in_str)
            event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
            desired=[3,4]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
            desired=[7,8,9]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
            x, y = map(event_cor[1][0],event_cor[0][0])
            ax = plt.gca()
            np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
            beach1 = beach(np1, xy=(x, y), width=0.09)
            ax.add_collection(beach1)
            for argv in sys.argv:
                if argv == "--topography":
                    try:
                        xpixels = 1000
                        map.arcgisimage(service='World_Shaded_Relief',
                                        xpixels = xpixels,
                                        verbose= False)
                    except:
                        pass

            plt.show()

            try:
                pathlist = Path(rel).glob('1-'+ str(sys.argv[5])+'.ASC')
            except:
                pathlist = Path(rel).glob('1-*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()
            map, x, y = make_map(data)

            mins = np.max(data[:,2])

            triang = tri.Triangulation(x, y)
            isbad = np.less(data_int, 0.01)
            mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang.set_mask(mask)
            plt.tricontourf(triang, data_int, cmap='YlOrRd')
            event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
            desired=[3,4]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
            desired=[7,8,9]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
            x, y = map(event_cor[1][0],event_cor[0][0])
            ax = plt.gca()
            np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
            beach1 = beach(np1, xy=(x, y), width=0.09)
            ax.add_collection(beach1)
            plt.colorbar()
            plt.title(path_in_str)
            for argv in sys.argv:
                if argv == "--topography":
                    try:
                        xpixels = 1000
                        map.arcgisimage(service='World_Shaded_Relief',
                                        xpixels = xpixels,
                                        verbose= False)
                    except:
                        pass

            plt.show()


def plot_integrated_timestep():

    evpath = 'events/' + str(sys.argv[1])
    C  = config.Config (evpath)
    Config = C.parseConfig ('config')
    cfg = ConfigObj (dict=Config)
    step = cfg.UInt ('step')
    step2 = cfg.UInt ('step_f2')

    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
            try:
                pathlist = Path(rel).glob('1-'+ str(sys.argv[5])+'.ASC')
            except:
                pathlist = Path(rel).glob('1-*boot*.ASC')
            maxs = 0.
            times = []
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    time = path_in_str[:-6]
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]
            try:
                pathlist = Path(rel).glob('1-'+ str(sys.argv[5])+'.ASC')
            except:
                pathlist = Path(rel).glob('1-*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):

                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()

            map, x, y = make_map(data)

            mins = np.max(data[:,2])

            triang = tri.Triangulation(x, y)
            isbad = np.less(data_int, 0.01)
            mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)

            triang.set_mask(mask)
            plt.tricontourf(triang, data_int, cmap='YlOrRd')

            plt.colorbar()
            plt.title(path_in_str)
            for argv in sys.argv:
                if argv == "--topography":
                    try:
                        xpixels = 1000
                        map.arcgisimage(service='World_Shaded_Relief',
                                        xpixels = xpixels,
                                        verbose= False)
                    except:
                        pass

            plt.show()

            pathlist = Path(rel).glob('0-*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()
            map, x, y = make_map(data)

            mins = np.max(data[:,2])
            import matplotlib.colors as colors
            import matplotlib.tri as tri

            triang = tri.Triangulation(x, y)
            isbad = np.less(data_int, 0.01)
            mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang.set_mask(mask)
            plt.tricontourf(triang, data_int, cmap='YlOrRd')

            plt.colorbar()
            plt.title(path_in_str)

            for argv in sys.argv:
                if argv == "--topography":
                    try:
                        xpixels = 1000
                        map.arcgisimage(service='World_Shaded_Relief',
                                        xpixels = xpixels,
                                        verbose= False)
                    except:
                        pass

            plt.show()

def plot_integrated_kite():
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            from kite import Scene
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
            pathlist = Path(rel).glob('0-9*.ASC')
            maxs = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]
            pathlist = Path(rel).glob('0-9*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()
            scd = Scene.load('/media/asteinbe/decepticon/playground/data/events/qaidam2009/insarnew/T319_20090708-20091021')

            data_dsc= scd.displacement
            eastings1 = np.arange(scd.frame.llLon,scd.frame.llLon+scd.frame.dE*scd.frame.cols,scd.frame.dE)
            northings1 = np.arange(scd.frame.llLat,scd.frame.llLat+scd.frame.dN*scd.frame.rows,scd.frame.dN)
            map = Basemap(projection='merc', llcrnrlon=num.min(eastings1),llcrnrlat=num.min(northings1),urcrnrlon=num.max(eastings1),urcrnrlat=num.max(northings1),
                          resolution='h', epsg=3395)
            ratio_lat = num.max(northings)/num.min(northings)
            ratio_lon = num.max(eastings)/num.min(eastings)

            map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)
            parallels = np.arange(num.min(northings1),num.max(northings),0.2)
            meridians = np.arange(num.min(eastings1),num.max(eastings),0.2)
            xpixels = 1000
           # map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            plt.tricontourf(x,y, data_int, cmap='hot', alpha=0.6)
            plt.colorbar()
            plt.title(path_in_str)
            ax = plt.gca()
            np1 = [101, 60, 83]
            x, y = map(95.76,37.64)

            beach1 = beach(np1, xy=(x, y), width=0.05)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            map.imshow(data_dsc)

            plt.show()

            pathlist = Path(rel).glob('1-9*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    if int(path_in_str[-7])==0:
                        data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                        data_int += np.nan_to_num(data[:,2])
                    else:
                        pass

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()
            map = Basemap(projection='merc', llcrnrlon=num.min(eastings1),llcrnrlat=num.min(northings1),urcrnrlon=num.max(eastings1),urcrnrlat=num.max(northings1),
                          resolution='h', epsg=3395)
            ratio_lat = num.max(northings)/num.min(northings)
            ratio_lon = num.max(eastings)/num.min(eastings)

            map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)
            parallels = np.arange(num.min(northings1),num.max(northings1),0.2)
            meridians = np.arange(num.min(eastings1),num.max(eastings1),0.2)
            xpixels = 1000
            #map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            plt.tricontourf(x,y, data_int, cmap='hot', alpha=0.6)
            plt.colorbar()
            plt.title(path_in_str)
            map.imshow(data_dsc)
            plt.show()

            pathlist = Path(rel).glob('1-9*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    if int(path_in_str[-7])>0:
                        data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                        data_int += np.nan_to_num(data[:,2])
                    else:
                        pass

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()
            map = Basemap(projection='merc', llcrnrlon=num.min(eastings1),llcrnrlat=num.min(northings1),urcrnrlon=num.max(eastings1),urcrnrlat=num.max(northings1),
                          resolution='h', epsg=3395)
            ratio_lat = num.max(northings)/num.min(northings)
            ratio_lon = num.max(eastings)/num.min(eastings)

            map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)
            parallels = np.arange(num.min(northings1),num.max(northings1),0.2)
            meridians = np.arange(num.min(eastings1),num.max(eastings1),0.2)
            xpixels = 1000
            #map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            plt.tricontourf(x,y, data_int, cmap='hot', alpha=0.6)
            plt.colorbar()
            plt.title(path_in_str)
            map.imshow(data_dsc)
            plt.show()


def plot_moving():
    datas = []
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
            pathlist = Path(rel).glob('**/1-0.1_*.ASC')
            for path in sorted(pathlist):
                try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    eastings = data[:,1]
                    northings =  data[:,0]
                except:
                    pass
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
            pathlist = Path(rel).glob('**/1-0.1_*.ASC')

        else:
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/' + str(sys.argv[3])
            pathlist = Path(rel).glob('**/*.ASC')
            for path in sorted(pathlist):
                try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    eastings = data[:,1]
                    northings =  data[:,0]
                except:
                    pass
            pathlist = Path(rel).glob('**/*.ASC')

        map, x, y = make_map(data)

        for path in sorted(pathlist):
            try:
                path_in_str = str(path)
                data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                eastings = data[:,1]
                northings =  data[:,0]
                x, y = map(data[:,1], data[:,0])
                datas.append(data[:,2])
            except:
                pass
        #mins = np.max(data[:,3])
        data = num.zeros(num.shape(data[:,2]))
        datas = num.asarray(datas)
        for dp in datas:
            data = data+dp
        scat = plt.tricontourf(x,y, data, vmin=num.max(data)*0.88)
        plt.show()

def plot_sembmax():
    evpath = 'events/'+ str(sys.argv[1])
    C  = config.Config (evpath)
    Config = C.parseConfig ('config')
    cfg = ConfigObj (dict=Config)
    filters = cfg.String('filters')
    filters = int(filters)
    step = cfg.UInt ('step')
    step2 = cfg.UInt ('step_f2')
    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    data, data_int, data_boot, data_int_boot, path_in_str, maxs, datamax = load(0)
    map, x, y = make_map(data)
    data = num.loadtxt(rel+'sembmax_0_boot0_P.txt', delimiter=' ')
    eastings = data[:,2]
    northings =  data[:,1]

    xpixels = 1000


    try:
        if sys.argv[3] is not None:
            with open(sys.argv[3], 'r') as fin:
                reader=csv.reader(fin)
                event_cor=[[(s) for s in row] for i,row in enumerate(reader)]
                event_cor_x = event_cor[:]
                event_cor_x = np.asarray(event_cor)
                event_cor_x = event_cor_x[:,1]
                event_cor_y = event_cor[:]
                event_cor_y = np.asarray(event_cor)
                event_cor_y = event_cor_y[:,2]
                list_x = []
                list_y = []
                for item in event_cor_x:
                    list_x.append(float(item))
                for item in event_cor_y:
                    list_y.append(float(item))
            x, y = map(list_x,list_y)
            ps = map.scatter(x,y,marker='o',c='k', s=size, cmap='seismic')
    except:
        pass

    ratio_lat = num.max(northings)/num.min(northings)
    ratio_lon = num.max(eastings)/num.min(eastings)
    #map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)
    event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
    desired=[3,4]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
    desired=[7,8,9]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
    x, y = map(event_cor[1][0],event_cor[0][0])
    ax = plt.gca()
    np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
    beach1 = beach(np1, xy=(x, y), width=1000.3, alpha=0.4)
    ax.add_collection(beach1)
    X,Y = np.meshgrid(eastings, northings)

    eastings, northings = map(X, Y)
    x, y = map(data[:,2], data[:,1])
    size =(data[:,3]/np.max(data[:,3]))*3000
    l = num.linspace(0,len(data[:,2])*step,len(data[:,2]))

    ps = map.scatter(x,y,marker='o',c=l, s=300, cmap='seismic')
    for i in range(0,len(x)):
        #if data[i,3]> np.max(data[:,3])*0.05:
            plt.text(x[i],y[i],'%s' %i)
    xpixels = 1000
#    map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
#    parallels = num.arange(num.min(northings),num.max(northings),ratio_lat)
#    meridians = num.arange(num.min(eastings),num.max(eastings),ratio_lon)
    #map.drawmeridians(meridians,labels=[1,1,1,1],linewidth=0.5, fontsize=10, dashes=[1,5])
    #map.drawparallels(parallels,labels=[1,1,1,1],linewidth=0.5, fontsize=10, dashes=[1,5])
    cbar = map.colorbar(ps,location='bottom',pad="5%", label='Time [s]')
    plt.show()
    try:
        rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
        data = num.loadtxt(rel+'sembmax_1_boot0_P.txt', delimiter=' ')
        eastings = data[:,2]
        northings =  data[:,1]

        map = Basemap(projection='merc', llcrnrlon=num.min(eastings),
                      llcrnrlat=num.min(northings),
                      urcrnrlon=num.max(eastings),
                      urcrnrlat=num.max(northings),
                      resolution='h', epsg=3395)

        try:
            if sys.argv[3] is not None:
                with open(sys.argv[3], 'r') as fin:
                    reader=csv.reader(fin)
                    event_cor=[[(s) for s in row] for i,row in enumerate(reader)]
                    event_cor_x = event_cor[:]
                    event_cor_x = np.asarray(event_cor)
                    event_cor_x = event_cor_x[:,1]
                    event_cor_y = event_cor[:]
                    event_cor_y = np.asarray(event_cor)
                    event_cor_y = event_cor_y[:,2]
                    list_x = []
                    list_y = []
                    for item in event_cor_x:
                        list_x.append(float(item))
                    for item in event_cor_y:
                        list_y.append(float(item))
                x, y = map(list_x,list_y)
                ps = map.scatter(x,y,marker='o',c='k', s=size, cmap='seismic')
        except:
            pass

        ratio_lat = num.max(northings)/num.min(northings)
        ratio_lon = num.max(eastings)/num.min(eastings)

        map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)

        event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
        desired=[3,4]
        with open(event, 'r') as fin:
            reader=csv.reader(fin)
            event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
        desired=[7,8,9]
        with open(event, 'r') as fin:
            reader=csv.reader(fin)
            event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
        x, y = map(event_cor[1][0],event_cor[0][0])
        ax = plt.gca()
        np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
        beach1 = beach(np1, xy=(x, y), width=1000.03, alpha=0.4)
        ax.add_collection(beach1)
        X,Y = np.meshgrid(eastings, northings)

        eastings, northings = map(X, Y)
        map.drawcoastlines(color='b',linewidth=1)

        x, y = map(data[:,2], data[:,1])
        for i in range(0,len(x)):
            if data[i,3]> np.max(data[:,3])*0.05:
                plt.text(x[i],y[i],'%s' %i)
        l = num.linspace(0,len(data[:,2])*step2,len(data[:,2]))
        size =(data[:,3]/np.max(data[:,3]))*3000
        ps = map.scatter(x,y,marker='o',c=l, s=size, cmap='seismic')
        xpixels = 1000
        map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
        parallels = num.arange(num.min(northings),num.max(northings),0.2)
        meridians = num.arange(num.min(eastings),num.max(eastings),0.2)
        cbar = map.colorbar(ps,location='bottom',pad="5%", label='Time [s]')
        plt.savefig(rel+'semblance_max_1.pdf', bbox_inches='tight')
        plt.show()
    except:
        pass


def plot_movingsembmax():
    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    data = num.loadtxt(rel+'sembmax_0_boot0_P.txt', delimiter=' ')
    eastings = data[:,2]
    northings =  data[:,1]
    xpixels = 1000
    map = Basemap(projection='merc', llcrnrlon=num.min(eastings),
                  llcrnrlat=num.min(northings),
                  urcrnrlon=num.max(eastings),
                  urcrnrlat=num.max(northings),
                  resolution='h', epsg=3395)
    ratio_lat = num.max(northings)/num.min(northings)
    ratio_lon = num.max(eastings)/num.min(eastings)

    map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)

    X,Y = np.meshgrid(eastings, northings)
    event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
    desired=[3,4]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
    desired=[7,8,9]
    with open(event, 'r') as fin:
        reader=csv.reader(fin)
        event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
    x, y = map(event_cor[1][0],event_cor[0][0])
    ax = plt.gca()
    np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
    beach1 = beach(np1, xy=(x, y), width=0.03, alpha=0.4)
    ax.add_collection(beach1)
    eastings, northings = map(X, Y)
    map.drawcoastlines(color='b',linewidth=1)
    map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
    parallels = num.arange(num.min(northings),num.max(northings),0.2)
    meridians = num.arange(num.min(eastings),num.max(eastings),0.2)
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
    map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
    x, y = map(data[:,2], data[:,1])
    size = num.shape(data[:,2])[0]
    l = range(0,size)
    si =(data[:,3]/np.max(data[:,3]))*300

    scat = map.scatter(x,y,marker='o',c=l, cmap='jet', s=si)
    axcolor = 'lightgoldenrodyellow'
    axamp = axes([0.2, 0.01, 0.65, 0.03])

    scorr = Slider(axamp, 'corr', 0, size, valinit=1)
    color=cm.rainbow(np.linspace(0,np.max(data[1,3]*1000),size))
    def update(val):
        corr = scorr.val
        i = int(corr)
        xx = np.vstack((x, y))
        scat.set_offsets(xx.T[i])
        scat.set_facecolor(color[int(data[i,3]*1000)])

        draw()

    scorr.on_changed(update)

    show(scat)
    try:
        rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
        data = num.loadtxt(rel+'sembmax_1_boot0_P.txt', delimiter=' ')
        eastings = data[:,2]
        northings =  data[:,1]
        xpixels = 1000
        map = Basemap(projection='merc', llcrnrlon=num.min(eastings),
                      llcrnrlat=num.min(northings),
                      urcrnrlon=num.max(eastings),
                      urcrnrlat=num.max(northings),
                      resolution='h', epsg=3395)
        ratio_lat = num.max(northings)/num.min(northings)
        ratio_lon = num.max(eastings)/num.min(eastings)

        map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)
        event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
        desired=[3,4]
        with open(event, 'r') as fin:
            reader=csv.reader(fin)
            event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
        desired=[7,8,9]
        with open(event, 'r') as fin:
            reader=csv.reader(fin)
            event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
        x, y = map(event_cor[1][0],event_cor[0][0])
        ax = plt.gca()
        np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
        beach1 = beach(np1, xy=(x, y), width=0.03, alpha=0.4)
        ax.add_collection(beach1)
        X,Y = np.meshgrid(eastings, northings)

        eastings, northings = map(X, Y)
        map.drawcoastlines(color='b',linewidth=1)
        map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
        parallels = num.arange(num.min(northings),num.max(northings),0.2)
        meridians = num.arange(num.min(eastings),num.max(eastings),0.2)
        map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
        map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
        x, y = map(data[:,2], data[:,1])
        size = num.shape(data[:,2])[0]
        l = range(0,size)
        si =(data[:,3]/np.max(data[:,3]))*300
        scat = map.scatter(x,y,marker='o',c=l, cmap='jet', s=si)
        axcolor = 'lightgoldenrodyellow'
        axamp = axes([0.2, 0.01, 0.65, 0.03])

        scorr = Slider(axamp, 'corr', 0, size, valinit=1)
        color=cm.rainbow(np.linspace(0,np.max(data[1,3]*1000),size))
        def update(val):
            corr = scorr.val
            i = int(corr)
            xx = np.vstack((x, y))
            scat.set_offsets(xx.T[i])
            scat.set_facecolor(color[int(data[i,3]*1000)])

            draw()

        scorr.on_changed(update)

        show(scat)
    except:
        pass

def sta_lta(data, dt, min_period):

    from scipy.signal import lfilter
    """
    The same STA/LTA as used in Flexwin.

    :copyright:
        Lion Krischer (krischer@geophysik.uni-muenchen.de), 2014
    :license:
        GNU General Public License, Version 3
        (http://www.gnu.org/copyleft/gpl.html)

    STA/LTA as used in FLEXWIN.

    :param data: The data array.
    :param dt: The sample interval of the data.
    :param min_period: The minimum period of the data.
    """
    Cs = 10 ** (-dt / min_period)
    Cl = 10 ** (-dt / (12 * min_period))
    TOL = 1e-9

    noise = data.max() / 1E5

    # 1000 samples should be more then enough to "warm up" the STA/LTA.
    extended_syn = np.zeros(len(data) + 1000, dtype=np.float64)
    # copy the original synthetic into the extended array, right justified
    # and add the noise level.
    extended_syn += noise
    extended_syn[-len(data):] += data

    # This piece of codes "abuses" SciPy a bit by "constructing" an IIR
    # filter that does the same as the decaying sum and thus avoids the need to
    # write the loop in Python. The result is a speedup of up to 2 orders of
    # magnitude in common cases without needing to write the loop in C which
    # would have a big impact in the ease of installation of this package.
    # Other than that its quite a cool little trick.
    a = [1.0, -Cs]
    b = [1.0]
    sta = lfilter(b, a, extended_syn)

    a = [1.0, -Cl]
    b = [1.0]
    lta = lfilter(b, a, extended_syn)

    # STA is now STA_LTA
    sta /= lta

    # Apply threshold to avoid division by very small values.
    sta[lta < TOL] = noise
    return sta[-len(data):]


def plot_semb():
    import matplotlib
    evpath = 'events/'+ str(sys.argv[1])
    C = config.Config(evpath)
    Config = C.parseConfig('config')
    cfg = ConfigObj(dict=Config)
    step = cfg.UInt('step')
    winlen = cfg.UInt('winlen')
    step2 = cfg.UInt('step_f2')
    winlen2 = cfg.UInt('winlen_f2')
    matplotlib.rcParams.update({'font.size': 22})
    rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
    astf = num.loadtxt(rel+'sembmax_0_boot0_P.txt', delimiter=' ')
    astf_data = astf[:, 3]
    trigger = sta_lta(astf_data, step, winlen)
    from scipy.signal import argrelextrema

    trigger[trigger<num.max(trigger*0.1)] =0
    extremas = argrelextrema(trigger, num.greater, order=4)
    minimas = argrelextrema(trigger, num.less, order=2)
    absmax = num.where(trigger>num.max(trigger)*0.2)
    print('duration from lf:', absmax[0][-1]*step-absmax[0][0]*step)

    fig = plt.figure()
    l = num.linspace(0,len(astf_data)*step,len(astf_data))
    plt.plot(l, astf_data ,'k')
    plt.ylabel('Semblance', fontsize=22)
    plt.xlabel('Time [s]', fontsize=22)

    plt.axvline(x=absmax[0][-1]*step, lw=4, c='r')
    plt.axvline(x=absmax[0][0]*step, lw=4, c='r')
    for ex in extremas[0]:
        plt.axvline(x=ex*step, lw=4, c='b')
    for ex in minimas[0]:
        plt.axvline(x=ex*step, lw=4, c='k')

    plt.savefig(rel+'semblance_0.pdf', bbox_inches='tight')
    plt.show()

    rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
    astf = num.loadtxt(rel+'sembmax_1_boot0_P.txt', delimiter=' ')
    astf_data = astf[:, 3]
    trigger = sta_lta(astf_data, step2, winlen2)
    trigger[trigger<num.max(trigger*0.1)] =0
    extremas = argrelextrema(trigger, num.greater, order=4)
    minimas = argrelextrema(trigger, num.less, order=2)
    absmax = num.where(trigger>num.max(trigger)*0.2)
    print('duration from hf:',absmax[0][-1]*step2-absmax[0][0]*step2)

    fig = plt.figure()
    l = num.linspace(0,len(astf_data)*step2,len(astf_data))
    plt.plot(l, astf_data ,'k')
    plt.ylabel('Semblance', fontsize=22)
    plt.xlabel('Time [s]', fontsize=22)

    plt.axvline(x=absmax[0][-1]*step2, lw=4, c='r')
    plt.axvline(x=absmax[0][0]*step2, lw=4, c='r')
    for ex in extremas[0]:
        plt.axvline(x=ex*step2, lw=4, c='b')
    for ex in minimas[0]:
        plt.axvline(x=ex*step2, lw=4, c='k')


    plt.savefig(rel+'semblance_1.pdf', bbox_inches='tight')
    plt.show()


def blobify():
    if len(sys.argv)<3:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':

            evpath = 'events/'+ str(sys.argv[1])
            C  = config.Config (evpath)
            Config = C.parseConfig ('config')
            cfg = ConfigObj (dict=Config)
            dimx = cfg.UInt ('dimx')
            dimy = cfg.UInt ('dimy')

            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

            import matplotlib.patches as mpatches

            from skimage.filters import threshold_otsu
            from skimage.segmentation import clear_border
            from skimage.measure import label, regionprops
            from skimage.morphology import closing, square
            from skimage.color import label2rgb

            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

            data_time = num.loadtxt(rel+'times_max_1_15.0.ASC', delimiter=' ')

            image = data_time[:,2]
            image = num.reshape(image, (dimx, dimy))
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.imshow(image)
            plt.show()

            # apply threshold
            thresh = threshold_otsu(image)
            bw = closing(image > thresh, square(3))

            # remove artifacts connected to image border
            cleared = clear_border(bw)

            # label image regions
            label_image = label(cleared)
            image_label_overlay = label2rgb(label_image, image=image)

            fig, ax = plt.subplots(figsize=(10, 6))
            ax.imshow(image_label_overlay)
            polygons = []
            for region in regionprops(label_image):
                # take regions with large enough areas
                if region.area >= 100:
                    # draw rectangle around segmented coins
                    minr, minc, maxr, maxc = region.bbox
                    rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,
                                              fill=False, edgecolor='red', linewidth=2)
                    ax.add_patch(rect)
                    polygons.append(Polygon((minc,maxc), (maxr,minr)))
            ax.set_axis_off()
            plt.tight_layout()
            plt.show()

            data = num.loadtxt(rel+'semb_cum_0_8.7.ASC', delimiter=' ')

            image = data[:,2]
            image = num.reshape(image, (dimx, dimy))
            # apply threshold
            thresh = threshold_otsu(image)
            bw = closing(image > thresh, square(3))

            # remove artifacts connected to image border
            cleared = clear_border(bw)

            # label image regions
            label_image = label(cleared)
            image_label_overlay = label2rgb(label_image, image=image)

            fig, ax = plt.subplots(figsize=(10, 6))
            ax.imshow(image_label_overlay)
            polygons = []
            for region in regionprops(label_image):
                # take regions with large enough areas
                if region.area >= 100:
                    # draw rectangle around segmented coins
                    minr, minc, maxr, maxc = region.bbox
                    rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,
                                              fill=False, edgecolor='red', linewidth=2)
                    ax.add_patch(rect)
                    polygons.append(Polygon((minc,maxc), (maxr,minr)))
            ax.set_axis_off()
            plt.tight_layout()
            plt.show()

            #load time for labeling
            labels = []
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.imshow(image_label_overlay)
            for polygon in polygons:
                l = 0
                for pe, pn in zip(data_time[:,0], data_time[:,1]):
                    point = Point(pn, pe)

                    if polygon.contains(point) is True:
                        labels.append(data_time[l,2])
                    l =+1
                    polygon.label= num.mean(labels)


def plot_scatter():
    if len(sys.argv)<3:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

            pathlist = Path(rel).glob('0-*.ASC')
            maxs = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    counter =+ 1
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]

            pathlist = Path(rel).glob('0-'+str(sys.argv[4])+('*.ASC'))
            data_int = num.zeros(num.shape(data[:, 2]))
            data_old = num.zeros(num.shape(data[:, 2]))
            time_grid = num.zeros(num.shape(data[:, 2]))
            counter = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])
                    for i in range(0, num.shape(data[:, 2])[0]):
                        if data[i,2] >= data_old[i]:
                            time_grid[i] = time_grid[i]+counter
                    data_old = np.nan_to_num(data[:,2])
                    counter =+ 1

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()

            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),
                          llcrnrlat=num.min(northings),
                          urcrnrlon=num.max(eastings),
                          urcrnrlat=num.max(northings),
                          resolution='h', epsg=3395)
            ratio_lat = num.max(northings)/num.min(northings)
            ratio_lon = num.max(eastings)/num.min(eastings)

            map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)

            xpixels = 1000

            eastings, northings = map(eastings, northings)
            parallels = num.arange(num.min(northings),num.max(northings),0.2)
            meridians = num.arange(num.min(eastings),num.max(eastings),0.2)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            data_old = num.zeros(num.shape(data[:, 2]))
            pathlist = Path(rel).glob('0-'+str(sys.argv[4])+('*.ASC'))
            i=0
            for path in sorted(pathlist):
                path_in_str = str(path)
                i = i+1

            colors = iter(cm.rainbow(np.linspace(0, 1, i)))
            cs =[]
            pathlist = Path(rel).glob('0-'+str(sys.argv[4])+('*.ASC'))
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    size =(data[:,2])*300
                    c = next(colors)
                    ps = map.scatter(x,y,marker='o',c=c, s=size, cmap='autumn_r')
                    data_old = np.nan_to_num(data[:,2])
                    cs.append(c)
            cmap_name = 'my_list'
            cms = LinearSegmentedColormap.from_list(
                cmap_name, cs, N=i)
            sm = plt.cm.ScalarMappable(cmap=cms, norm=plt.Normalize(vmin=0, vmax=1))
            sm._A = []
            cb = plt.colorbar(sm, orientation="horizontal")
            cb.outline.set_visible(False)
            cb.set_ticks([])
            cb.set_label('Time ->',fontsize=22)
            plt.title(path_in_str)
            ax = plt.gca()
            event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
            desired=[3,4]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
            desired=[7,8,9]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
            x, y = map(event_cor[1][0],event_cor[0][0])
            ax = plt.gca()
            np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
            beach1 = beach(np1, xy=(x, y), width=0.03, alpha=0.4)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

            plt.show()

            pathlist = Path(rel).glob('1-*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()


            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),
                          llcrnrlat=num.min(northings),
                          urcrnrlon=num.max(eastings),
                          urcrnrlat=num.max(northings),
                          resolution='h', epsg=3395)
            ratio_lat = num.max(northings)/num.min(northings)
            ratio_lon = num.max(eastings)/num.min(eastings)

            map.drawmapscale(num.min(eastings)+ratio_lon*0.25, num.min(northings)+ratio_lat*0.25, num.mean(eastings), num.mean(northings), 30)

            xpixels = 1000

            eastings, northings = map(eastings, northings)
            parallels = num.arange(num.min(northings),num.max(northings),0.2)
            meridians = num.arange(num.min(eastings),num.max(eastings),0.2)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            data_old = num.zeros(num.shape(data[:, 2]))
            pathlist = Path(rel).glob('1-'+str(sys.argv[4])+('*.ASC'))
            i=0
            for path in sorted(pathlist):
                path_in_str = str(path)
                i = i+1

            colors = iter(cm.rainbow(np.linspace(0, 1, i)))
            cs =[]
            pathlist = Path(rel).glob('1-'+str(sys.argv[4])+('*.ASC'))
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    size =(data[:,2])*300
                    c = next(colors)
                    ps = map.scatter(x,y,marker='o',c=c, s=size, cmap='autumn_r')
                    data_old = np.nan_to_num(data[:,2])
                    cs.append(c)
            cmap_name = 'my_list'
            cms = LinearSegmentedColormap.from_list(
                cmap_name, cs, N=i)
            sm = plt.cm.ScalarMappable(cmap=cms, norm=plt.Normalize(vmin=0, vmax=1))
            sm._A = []
            cb = plt.colorbar(sm, orientation="horizontal")
            cb.outline.set_visible(False)
            cb.set_ticks([])
            cb.set_label('Time ->',fontsize=22)
            plt.title(path_in_str)
            ax = plt.gca()
            event = 'events/'+ str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
            desired=[3,4]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
            desired=[7,8,9]
            with open(event, 'r') as fin:
                reader=csv.reader(fin)
                event_mech=[[float(s[-3:]) for s in row] for i,row in enumerate(reader) if i in desired]
            x, y = map(event_cor[1][0],event_cor[0][0])
            ax = plt.gca()
            np1 = [event_mech[0][0], event_mech[1][0], event_mech[2][0]]
            beach1 = beach(np1, xy=(x, y), width=0.03, alpha=0.4)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

            plt.show()


def empiricial_timeshifts():
        import _pickle as pickle

        evpath = 'events/'+ str(sys.argv[1])
        C = config.Config(evpath)
        Config = C.parseConfig('config')
        cfg = ConfigObj(dict=Config)
        sembpath = evpath + '/work/semblance'
        stations = []
        refs = []
        rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

        if cfg.Bool('synthetic_test') is True:
            Syn_in = C.parseConfig('syn')
            syn_in = SynthCfg(Syn_in)
            lat_ev = float(syn_in.lat_0())
            lon_ev = float(syn_in.lon_0())
        else:
            event = 'events/' + str(sys.argv[1]) + '/' + str(sys.argv[1])+'.origin'
            desired=[3,4]
            with open(event, 'r') as fin:
                reader = csv.reader(fin)
                event_cor=[[float(s[6:]) for s in row] for i,row in enumerate(reader) if i in desired]
                lat_ev, lon_ev = event_cor[1][0], event_cor[0][0]


        pathlist = Path(rel).glob('*.shift*')
        for path in sorted(pathlist):
                path_in_str = str(path)
                if path_in_str[-1] != "s":
                    f = open(path_in_str, 'rb')
                    refshifts = pickle.load(f)
                    f.close()
                    for s in refshifts.values():
                        refs.append(s)

                else:
                    f = open(path_in_str, 'rb')
                    refshifts_stations = pickle.load(f)
                    f.close()
                    for s in refshifts_stations.values():
                        stations.append(s)
        bazis = []
        dists = []
        for s in stations:
            print(s)
            b = orthodrome.azimuth(s[0], s[1], lat_ev, lon_ev)
            dists.append(orthodrome.distance_accurate50m(s[0], s[1], lat_ev, lon_ev))

            if b>=0.:
                bazi = b
            elif b<0.:
                bazi = 360.+b
            bazis.append(bazi)
        plt.figure()
        plt.scatter(refs, bazis)
        plt.show()
        plt.figure()
        plt.scatter(refs, dists)
        plt.show()

def main():
    if len(sys.argv)<3:
        print("input: eventname plot_name,\
         available plot_name: movie, sembmax, nce, interactive_max, cluster")
    else:
        event = sys.argv[1]
        if sys.argv[2] == 'movie':
            plot_movie()
        elif sys.argv[2] == 'sembmax':
            plot_sembmax()
        elif sys.argv[2] == 'semblance':
            plot_semb()
        elif sys.argv[2] == 'interactive_max':
            plot_movingsembmax()
        elif sys.argv[2] == 'cluster':
            plot_cluster()
        elif sys.argv[2] == 'moving':
            plot_moving()
        elif sys.argv[2] == 'integrated':
            plot_integrated()
        elif sys.argv[2] == 'integrated_kite':
            plot_integrated_kite()
        elif sys.argv[2] == 'integrated_scatter':
            integrated_scatter()
        elif sys.argv[2] == 'scatter':
            plot_scatter()
        elif sys.argv[2] == 'beampower':
            beampower()
        elif sys.argv[2] == 'blobify':
            blobify()
        elif sys.argv[2] == 'inspect_spectrum':
            inspect_spectrum()
        elif sys.argv[2] == 'semb_equal':
            plot_semb_equal()
        elif sys.argv[2] == 'integrated_timestep':
            plot_integrated_timestep()
        elif sys.argv[2] == 'distance_time':
            distance_time()
        elif sys.argv[2] == 'integrated_movie':
            plot_integrated_movie()
        elif sys.argv[2] == 'timeshifts':
            empiricial_timeshifts()
        elif sys.argv[2] == 'distance_time_bootstrap':
            distance_time_bootstrap()
        elif sys.argv[2] == 'timeshifts_map':
            plot_timeshift_map()
