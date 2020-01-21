from pyrocko.model import event
from pyrocko import orthodrome, util
from pyrocko.guts import dump
from pyrocko.model import Geometry
from palantiri.common import ConfigFile
from palantiri.common.ConfigFile import ConfigObj, FilterCfg, OriginCfg, SynthCfg
from palantiri.common import Globals
from palantiri.tools import config
import sys
from pathlib import Path
import numpy as num
global evpath

def duplicate_property(array):
    ndims = len(array.shape)
    if ndims == 1:
        return num.hstack((array, array))
    elif ndims == 2:
        return num.vstack((array, array))
    else:
        raise TypeError('Only 1-2d data supported!')



def load(filter, step=None, path=None):
            if path is None:
                rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
            else:
                rel = path
            boot = False
            if path is not None:
                evpath = path
            else:
                evpath = 'events/'+ str(sys.argv[1])
            C  = config.Config (evpath)
            Config = C.parseConfig ('config')
            cfg = ConfigObj (dict=Config)
            dimx = int(Config['dimx'])
            dimy = int(Config['dimy'])
            data_int = None
            data = None
            data_boot = None
            data_int_boot = None
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
                    pathlist = Path(rel).glob('*.ASC')
                except:
                    pathlist = Path(rel).glob('%s-*%s.ASC' % (filter,phase))
            else:
                try:
                    try:
                        pathlist = Path(rel).glob('*0%s.ASC' % step)
                    except:
                        pathlist = Path(rel).glob('*%s.ASC' % step)
                except:
                    pathlist = Path(rel).glob('%s-*00%s_*%s.ASC' % (filter, step, phase))
            maxs = 0.
            count = 0
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    maxd = num.max(data[:, 2])
                    count = count + 1
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
                        for k in num.nan_to_num(data[:,2]):
                            if k>data_int[i]:
                                data_int[i]= k
                            if num.max(datamax) == 0:
                                data_int[i]= 0
                            i = i+1
                try:
                    if sys.argv[4] == 'boot':
                        boot = True
                        if step is None:
                            try:
                                pathlist = Path(rel).glob('%s-*boot*'+ str(sys.argv[5])+'*.ASC' % filter)
                            except:
                                pathlist = Path(rel).glob('%s-*boot*%s.ASC' % (filter, phase))
                        else:
                            try:
                                pathlist = Path(rel).glob('%s-*boot*'+ str(sys.argv[5])+'00%s_*.ASC' % (filter, step))
                            except:
                                pathlist = Path(rel).glob('%s-*boot00%s_*%s.ASC' % (filter, step, phase))
                        data_int_boot = num.zeros(num.shape(data[:, 2]))
                        for path in sorted(pathlist):
                                path_in_str = str(path)
                                data_boot = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                                i = 0
                                for k in num.nan_to_num(data[:,2]):
                                    if k>data_int_boot[i]:
                                        data_int_boot[i]= k
                                    if num.max(datamax) == 0:
                                        data_int[i]= 0
                                    i = i+1
                except IndexError:
                    pass

            if sys.argv[3] == 'combined':
                if step is None:
                    try:
                        pathlist = Path(rel).glob('%s-'+ str(sys.argv[5])+'*.ASC' % filter)
                    except:
                        pathlist = Path(rel).glob('%s*-%s*.ASC' % (filter,phase))
                else:
                    try:
                        pathlist = Path(rel).glob('%s-'+ str(sys.argv[5])+'00%s_*.ASC' % (filter, step))
                    except:
                        pathlist = Path(rel).glob('%s-*00%s_*%s.ASC' % (filter, step, phase))
                data_int = num.zeros(num.shape(data[:, 2]))
                for path in sorted(pathlist):
                        path_in_str = str(path)
                        if path_in_str[-14] is not "o":
                            data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                            data_int += num.nan_to_num(data[:,2])

                try:
                    if sys.argv[4] == 'boot':
                        boot = True

                        if step is None:
                            try:
                                pathlist = Path(rel).glob('%s-*boot*'+ str(sys.argv[5])+'*.ASC' % filter)
                            except:
                                pathlist = Path(rel).glob('%s-*boot*.ASC' % filter)
                        else:
                            try:
                                pathlist = Path(rel).glob('%s-*boot*'+ str(sys.argv[5])+'00%s_*.ASC' % (filter, step))
                            except:
                                pathlist = Path(rel).glob('%s-*boot*00%s_*.ASC' % (filter, step))
                        data_int_boot = num.zeros(num.shape(data[:, 2]))
                        for path in sorted(pathlist):
                                path_in_str = str(path)
                                data_boot = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                                data_int_boot += num.nan_to_num(data_boot[:,2])
                except IndexError:
                    pass
            return data, data_int, data_boot, data_int_boot, path_in_str, maxs, datamax, count



def from_palantiri():
    km = 1000.
    try:
        path = sys.argv[3]
        evpath = path
    except:
        path = None
        evpath = 'events/'+ str(sys.argv[1])


    C  = config.Config(evpath)
    Origin = C.parseConfig('origin')
    Config = C.parseConfig('config')
    cfg = ConfigObj(dict=Config)
    step = cfg.UInt('step')
    step2 = cfg.UInt('step_f2')
    duration = cfg.UInt('duration')
    forerun = cfg.UInt('forerun')
    deltat = step
    deltat2 = step2
    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

    dimx = int(Config['dimx'])
    dimy = int(Config['dimy'])

    origin = OriginCfg(Origin)
    depth = origin.depth()*1000.
    ev = event.Event(lat=origin.lat(), lon=origin.lon(), depth=depth, time=util.str_to_time(origin.time()))
    data, data_int, data_boot, data_int_boot, path_in_str, maxs, datamax, n_files = load(0, path=path)
    values_orig = data[:, 2]
    #values_orig = num.append(values_orig, num.array([0., 0.]))

    lat_orig = data[:, 1]
    lon_orig = data[:, 0]

    ncorners = 4
    lon_grid_orig = num.linspace(num.min(lat_orig), num.max(lat_orig), (dimy))
    lat_grid_orig = num.linspace(num.min(lon_orig), num.max(lon_orig), dimx)

    if path is None:
        ntimes = int((forerun+duration)/step)
    else:
        ntimes = n_files

    verts = []
    lon_diff = ((lon_orig)[dimy+1]-(lon_orig)[0])/4.
    lat_diff = ((lat_orig)[1]-(lat_orig)[0])/4.

    dist = orthodrome.distance_accurate50m(lat_grid_orig[1], lon_grid_orig[1], lat_grid_orig[0], lon_grid_orig[0])

    for x,y in zip(lon_orig, lat_orig):

            xyz = ([dist/2.8, dist/2.8, depth], [-dist/2.8, dist/2.8, depth],[-dist/2.8, -dist/2.8, depth], [dist/2.8, -dist/2.8, depth] )
            latlon = ([x,y], [x,y], [x,y], [x,y])
            patchverts = num.hstack((latlon, xyz))
            verts.append(patchverts)


    vertices = num.vstack(verts)

    npatches = int(len(vertices)) #*2?
    faces1 = num.arange(ncorners * npatches, dtype='int64').reshape(
        npatches, ncorners)
    faces2 = num.fliplr(faces1)
    faces = num.vstack((faces2, faces1))
    srf_semblance_list = []
    for i in range(0,ntimes):
        if len(sys.argv)<4:
            print("missing input arrayname")
        else:
                data, data_int, data_boot, data_int_boot, path_in_str, maxsb, datamaxb, n_files = load(0, step=i, path=path)
                srf_semblance = data[:,2]
                #srf_semblance = num.append(srf_semblance, num.array([0., 0.]))
                srf_semblance = duplicate_property(srf_semblance)
                srf_semblance_list.append(srf_semblance)
                print(len(srf_semblance))
    srf_semblance = num.asarray(srf_semblance_list).T
    srf_times = num.linspace(0, forerun+duration, ntimes)
    geom = Geometry(times=srf_times, event=ev)
    geom.setup(vertices, faces)
    sub_headers = tuple([str(i) for i in srf_times])
    geom.add_property((('semblance', 'float64', sub_headers)), srf_semblance)
    dump(geom, filename='geom.yaml')


def main():
    if len(sys.argv)<3:
        print("input: eventname plot_name,\
         available plot_name: movie, sembmax, nce, interactive_max, cluster")
    else:
        event = sys.argv[1]
    if sys.argv[2] == 'export_geometry':
        from_palantiri()
