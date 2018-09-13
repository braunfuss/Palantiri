from affine import Affine
from pyproj import Proj, transform
from mpl_toolkits.basemap import Basemap
import numpy as num
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
from pathlib import Path
import sys
from matplotlib.pyplot import cm
from matplotlib.widgets import Slider
from pylab import plot, show, figure, scatter, axes, draw
from itertools import cycle
import random
import csv
from obspy.imaging.beachball import beach

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
    map.drawcoastlines()
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
        data = num.loadtxt(path_in_str, delimiter=' ', usecols=(0,2,3))
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


    plt.show()

def plot_movie():
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
            pathlist = Path(rel).glob('*.ASC')
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
                    map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
                            resolution='h')
                    parallels = np.arange(num.min(northings),num.max(northings),0.2)
                    meridians = np.arange(num.min(eastings),num.max(eastings),0.2)
                    xpixels = 1000
                    #map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
                    eastings, northings = map(eastings, northings)
                    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
                    map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
                    x, y = map(data[::10,1], data[::10,0])
                    mins = np.max(data[:,2])
                    plt.tricontourf(x,y, data[::10,2], cmap='hot', vmin=0., vmax=maxs)
                    plt.colorbar()
                    plt.title(path_in_str)
                    plt.savefig(path_in_str+'.pdf', bbox_inches='tight')
                    plt.close()
            #    except:
                #    plt.close()
                #    pass

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
                    map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
                            resolution='h')
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
                    plt.savefig(path_in_str+'.pdf', bbox_inches='tight')
                    plt.close()
                except:
                    plt.close()
                    pass


def plot_scatter():
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work-a/semblance/'
            import matplotlib
            matplotlib.rcParams.update({'font.size': 32})
            pathlist = Path(rel).glob('1-13*.ASC')
            maxs = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]
            pathlist = Path(rel).glob('1-13*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()
	    print(np.shape(data[:,2]))
#            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
#                    resolution='h',epsg = 4269)
            map = Basemap( projection='cyl',\
                    llcrnrlon=95.25, \
                    llcrnrlat=37, \
                    urcrnrlon=97.25, \
                    urcrnrlat=38, \
                    resolution='h',epsg = 4269)
            parallels = np.arange(37,38,1.)
            meridians = np.arange(95.5,97.5,0.5)
            xpixels = 1000
           # map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])

            #l = range(0,num.shape(data[:,2])[0])
            l = [i for i in range(1600) for _ in range(70)]
            l = sorted(range(100)*16)
            size = (data[:,2]/np.max(data[:,2]))*300
            ps = map.scatter(x,y,marker='o',c=l, s=size, cmap='autumn_r')

        #    for i in range(0,len(x)):
        #        if data[i,2]> np.max(data[:,2])*0.05:
        #            plt.text(x[i],y[i],'%s' %i)
            #data_int[data_int<np.max(data_int)*0.000001]=np.nan

            #plt.tricontourf(x,y, data_int, cmap='hot',alpha=0.6)

            #plt.tricontourf(x,y, data_int, cmap='hot',norm=colors.Normalize(vmin=0.1, vmax=1.1))
            plt.colorbar(orientation="horizontal")
            plt.title(path_in_str)
            ax = plt.gca()
            np1 = [101, 60, 83]
            x, y = map(95.76,37.64)

            beach1 = beach(np1, xy=(x, y), width=0.05)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

            plt.show()

            pathlist = Path(rel).glob('0-6*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()

#            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
#                    resolution='h',epsg = 4269)
            map = Basemap( projection='cyl',\
                    llcrnrlon=95.25, \
                    llcrnrlat=37, \
                    urcrnrlon=97.25, \
                    urcrnrlat=38, \
                    resolution='h',epsg = 4269)
            parallels = np.arange(37,38,1.)
            meridians = np.arange(95.5,97.5,0.5)
            xpixels = 1000
           # map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])

            #l = range(0,num.shape(data[:,2])[0])
            l = [i for i in range(1600) for _ in range(70)]
            l = sorted(range(160)*10)
            size = (data[:,2]/np.max(data[:,2]))*300
            ps = map.scatter(x,y,marker='o',c=l, s=size, cmap='winter_r')

    #        for i in range(0,len(x)):
    #            if data[i,2]> np.max(data[:,2])*0.05:
    #                plt.text(x[i],y[i],'%s' %i)
            #data_int[data_int<np.max(data_int)*0.000001]=np.nan

            #plt.tricontourf(x,y, data_int, cmap='hot',alpha=0.6)

            #plt.tricontourf(x,y, data_int, cmap='hot',norm=colors.Normalize(vmin=0.1, vmax=1.1))
            plt.colorbar(orientation="horizontal")
            plt.title(path_in_str)
            ax = plt.gca()
            np1 = [101, 60, 83]
            x, y = map(95.76,37.64)

            beach1 = beach(np1, xy=(x, y), width=0.05)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

            plt.show()

def plot_integrated():
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

            pathlist = Path(rel).glob('0-16*.ASC')
            maxs = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]
            pathlist = Path(rel).glob('0-16*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            print np.shape(data_int)

            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    i = 0
                    for k in np.nan_to_num(data[:,2]):
                        if k>data_int[i]:
                            data_int[i]= k
                        i = i+1
                    #data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()

#            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
#                    resolution='h',epsg = 4269)
            map = Basemap( projection='cyl',\
                    llcrnrlon=95.25, \
                    llcrnrlat=37, \
                    urcrnrlon=97.25, \
                    urcrnrlat=38, \
                    resolution='h',epsg = 4269)
            parallels = np.arange(37,38,1.)
            meridians = np.arange(95.5,97.5,0.5)
            xpixels = 1000
           # map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            #data_int[data_int<np.max(data_int)*0.000001]=np.nan
            import matplotlib.colors as colors
            import matplotlib.tri as tri
            #mask = np.ma.masked_where(data_int < 0.4, data_int)
            #mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang = tri.Triangulation(x, y)
            isbad = np.less(data_int, 0.085)
            mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            levels = np.arange(0., 1.05, 0.025)
            triang.set_mask(mask)
            plt.tricontourf(triang, data_int, cmap='cool')
            #plt.tricontourf(x,y, data_int, cmap='hot',alpha=0.6)

            #plt.tricontourf(x,y, data_int, cmap='hot',norm=colors.Normalize(vmin=0.1, vmax=1.1))
            plt.colorbar(orientation="horizontal")
            plt.title(path_in_str)
            ax = plt.gca()
            np1 = [116, 61, 91]
            x, y = map(96.476,37.529)

            beach1 = beach(np1, xy=(x, y), width=0.05)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

            plt.show()

            pathlist = Path(rel).glob('0-13*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()

#            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
#                    resolution='h',epsg = 4269)
            map = Basemap( projection='cyl',\
                    llcrnrlon=95.25, \
                    llcrnrlat=37, \
                    urcrnrlon=97.25, \
                    urcrnrlat=38, \
                    resolution='h',epsg = 4269)
            parallels = np.arange(37,38,1.)
            meridians = np.arange(95.5,97.5,0.5)
            xpixels = 1000
           # map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            #data_int[data_int<np.max(data_int)*0.000001]=np.nan
            import matplotlib.colors as colors
            import matplotlib.tri as tri
            #mask = np.ma.masked_where(data_int < 0.4, data_int)
            #mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang = tri.Triangulation(x, y)
            isbad = np.less(data_int, 0.01)
            mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang.set_mask(mask)
            plt.tricontourf(triang, data_int, cmap='YlOrRd')
            #plt.tricontourf(x,y, data_int, cmap='hot',alpha=0.6)

            #plt.tricontourf(x,y, data_int, cmap='hot',norm=colors.Normalize(vmin=0.1, vmax=1.1))
            plt.colorbar()
            plt.title(path_in_str)
            ax = plt.gca()
            np1 = [116, 61, 91]
            x, y = map(96.476,37.529)

            beach1 = beach(np1, xy=(x, y), width=0.05)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

            plt.show()



def plot_integrated_timestep():
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

            pathlist = Path(rel).glob('1-13*.ASC')
            maxs = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]
            pathlist = Path(rel).glob('1-13*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()

#            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
#                    resolution='h',epsg = 4269)
            map = Basemap( projection='cyl',\
                    llcrnrlon=95.25, \
                    llcrnrlat=37, \
                    urcrnrlon=97.25, \
                    urcrnrlat=38, \
                    resolution='h',epsg = 4269)
            parallels = np.arange(37,38,1.)
            meridians = np.arange(95.5,97.5,0.5)
            xpixels = 1000
           # map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            #data_int[data_int<np.max(data_int)*0.000001]=np.nan
            import matplotlib.colors as colors
            import matplotlib.tri as tri
            #mask = np.ma.masked_where(data_int < 0.4, data_int)
            #mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang = tri.Triangulation(x, y)
            isbad = np.less(data_int, 0.01)
            mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)

            triang.set_mask(mask)
            plt.tricontourf(triang, data_int, cmap='YlOrRd')
            #plt.tricontourf(x,y, data_int, cmap='hot',alpha=0.6)

            #plt.tricontourf(x,y, data_int, cmap='hot',norm=colors.Normalize(vmin=0.1, vmax=1.1))
            plt.colorbar()
            plt.title(path_in_str)
            ax = plt.gca()
            np1 = [116, 61, 91]
            x, y = map(96.476,37.529)

            beach1 = beach(np1, xy=(x, y), width=0.05)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

            plt.show()

            pathlist = Path(rel).glob('0-13*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()

#            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
#                    resolution='h',epsg = 4269)
            map = Basemap( projection='cyl',\
                    llcrnrlon=95.25, \
                    llcrnrlat=37, \
                    urcrnrlon=97.25, \
                    urcrnrlat=38, \
                    resolution='h',epsg = 4269)
            parallels = np.arange(37,38,1.)
            meridians = np.arange(95.5,97.5,0.5)
            xpixels = 1000
           # map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            #data_int[data_int<np.max(data_int)*0.000001]=np.nan
            import matplotlib.colors as colors
            import matplotlib.tri as tri
            #mask = np.ma.masked_where(data_int < 0.4, data_int)
            #mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang = tri.Triangulation(x, y)
            isbad = np.less(data_int, 0.01)
            mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang.set_mask(mask)
            plt.tricontourf(triang, data_int, cmap='YlOrRd')
            #plt.tricontourf(x,y, data_int, cmap='hot',alpha=0.6)

            #plt.tricontourf(x,y, data_int, cmap='hot',norm=colors.Normalize(vmin=0.1, vmax=1.1))
            plt.colorbar()
            plt.title(path_in_str)
            ax = plt.gca()
            np1 = [116, 61, 91]
            x, y = map(96.476,37.529)

            beach1 = beach(np1, xy=(x, y), width=0.05)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

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
                    resolution='h',epsg = 4269)
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
                    resolution='h',epsg = 4269)
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
                    resolution='h',epsg = 4269)
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

        map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
                resolution='h')
        parallels = np.arange(num.min(northings),num.max(northings),0.2)
        meridians = np.arange(num.min(eastings),num.max(eastings),0.2)

        eastings, northings = map(eastings, northings)
        map.drawcoastlines(color='b',linewidth=3)
        map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
        map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
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
    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    data = num.loadtxt(rel+'sembmax_0.txt', delimiter=' ')
    eastings = data[:,2]
    northings =  data[:,1]

    map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
            resolution='h',epsg = 4269)
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

    x, y = map(data[:,2], data[:,1])
    l = range(0,num.shape(data[:,2])[0])
    size = (data[:,3]/np.max(data[:,3]))*300
    ps = map.scatter(x,y,marker='o',c=l, s=size, cmap='seismic')
    for i in range(0,len(x)):
        if data[i,3]> np.max(data[:,3])*0.05:
            plt.text(x[i],y[i],'%s' %i)
    xpixels = 1000
    map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
    parallels = num.arange(num.min(northings),num.max(northings),0.2)
    meridians = num.arange(num.min(eastings),num.max(eastings),0.2)
    map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
    map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
    cbar = map.colorbar(ps,location='bottom',pad="5%", label='Time [s]')
    plt.savefig(rel+'semblance_max_0.pdf', bbox_inches='tight')
    plt.show()
    try:
        rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
        data = num.loadtxt(rel+'sembmax_1.txt', delimiter=' ')
        eastings = data[:,2]
        northings =  data[:,1]

        map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
                resolution='h',epsg = 4269)

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

        x, y = map(data[:,2], data[:,1])
        for i in range(0,len(x)):
            if data[i,3]> np.max(data[:,3])*0.05:
                plt.text(x[i],y[i],'%s' %i)
        l = range(0,num.shape(data[:,2])[0])
        size = (data[:,3]/np.max(data[:,3]))*3000
        ps = map.scatter(x,y,marker='o',c=l, s=size, cmap='seismic')
        xpixels = 1000
        map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
        parallels = num.arange(num.min(northings),num.max(northings),0.2)
        meridians = num.arange(num.min(eastings),num.max(eastings),0.2)
        map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
        map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
        cbar = map.colorbar(ps,location='bottom',pad="5%", label='Time [s]')
        plt.savefig(rel+'semblance_max_1.pdf', bbox_inches='tight')
        plt.show()
    except:
        pass


def plot_movingsembmax():
    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    data = num.loadtxt(rel+'sembmax_0.txt', delimiter=' ')
    eastings = data[:,2]
    northings =  data[:,1]
    xpixels = 1000
    map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
            resolution='h',epsg = 4269)

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
    si = (data[:,3]/np.max(data[:,3]))*300

    scat = map.scatter(x,y,marker='o',c=l, cmap='jet', s=si)
    axcolor = 'lightgoldenrodyellow'
    axamp = axes([0.2, 0.01, 0.65, 0.03])

    scorr = Slider(axamp, 'corr', 0, size, valinit=1)
    color=cm.rainbow(np.linspace(0,np.max(data[1,3]*1000),size))
    def update(val):
        corr = scorr.val
        i = int(corr)
        xx = np.vstack ((x, y))
        scat.set_offsets (xx.T[i])
        scat.set_facecolor(color[int(data[i,3]*1000)])

        draw()

    scorr.on_changed(update)

    show(scat)
    try:
        rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
        data = num.loadtxt(rel+'sembmax_1.txt', delimiter=' ')
        eastings = data[:,2]
        northings =  data[:,1]
        xpixels = 1000
        map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
                resolution='h',epsg = 4269)
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
        si = (data[:,3]/np.max(data[:,3]))*300
        scat = map.scatter(x,y,marker='o',c=l, cmap='jet', s=si)
        axcolor = 'lightgoldenrodyellow'
        axamp = axes([0.2, 0.01, 0.65, 0.03])

        scorr = Slider(axamp, 'corr', 0, size, valinit=1)
        color=cm.rainbow(np.linspace(0,np.max(data[1,3]*1000),size))
        def update(val):
            corr = scorr.val
            i = int(corr)
            xx = np.vstack ((x, y))
            scat.set_offsets (xx.T[i])
            scat.set_facecolor(color[int(data[i,3]*1000)])

            draw()

        scorr.on_changed(update)

        show(scat)
    except:
        pass


def plot_semb():
    import matplotlib
    matplotlib.rcParams.update({'font.size': 22})
    rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
    astf = num.loadtxt(rel+'sembmax_0.txt', delimiter=' ')
    astf_data = astf[:, 3]
    fig = plt.figure()
    plt.plot(astf_data ,'k')
    plt.ylabel('Semblance', fontsize=22)
    plt.xlabel('Time [s]', fontsize=22)
    plt.savefig(rel+'semblance_0.pdf', bbox_inches='tight')
    plt.show()
    try:
        rel = 'events/' + str(sys.argv[1]) + '/work/semblance/'
        astf = num.loadtxt(rel+'sembmax_1.txt', delimiter=' ')
        fig = plt.figure()
        astf_data = astf[:, 3]

        plt.plot(astf_data, 'k')
        plt.ylabel('Beampower', fontsize=22)

        plt.xlabel('Time [s]', fontsize=22)

        plt.savefig(rel+'semblance_1.pdf', bbox_inches='tight')
        plt.show()
    except:
        pass


def integrated_scatter():
    if len(sys.argv)<4:
        print("missing input arrayname")
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'

            pathlist = Path(rel).glob('1-13*.ASC')
            maxs = 0.
            for path in sorted(pathlist):
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    max = np.max(data[:, 2])
                    if maxs < max:
                        maxs = max
                        datamax = data[:, 2]
            pathlist = Path(rel).glob('1-13*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()

#            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
#                    resolution='h',epsg = 4269)
            map = Basemap( projection='cyl',\
                    llcrnrlon=95.25, \
                    llcrnrlat=37, \
                    urcrnrlon=97.25, \
                    urcrnrlat=38, \
                    resolution='h',epsg = 4269)
            parallels = np.arange(37,38,1.)
            meridians = np.arange(95.5,97.5,0.5)
            xpixels = 1000
           # map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            #data_int[data_int<np.max(data_int)*0.000001]=np.nan
            #import matplotlib.colors as colors
            #import matplotlib.tri as tri
            #mask = np.ma.masked_where(data_int < 0.4, data_int)
            #mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)


            #plt.tricontourf(triang, data_int, cmap='YlOrRd')
		            #plt.tricontourf(x,y, data_int, cmap='hot',alpha=0.6)
	    plt.scatter(x, y, data[:,2])
            #plt.tricontourf(x,y, data_int, cmap='hot',norm=colors.Normalize(vmin=0.1, vmax=1.1))
            plt.colorbar()
            plt.title(path_in_str)
            ax = plt.gca()
            np1 = [116, 61, 91]
            x, y = map(96.476,37.529)

            beach1 = beach(np1, xy=(x, y), width=0.05)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

            plt.show()

            pathlist = Path(rel).glob('0-13*.ASC')
            data_int = num.zeros(num.shape(data[:, 2]))
            for path in sorted(pathlist):
            #    try:
                    path_in_str = str(path)
                    data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)
                    data_int += np.nan_to_num(data[:,2])

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()

#            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
#                    resolution='h',epsg = 4269)
            map = Basemap( projection='cyl',\
                    llcrnrlon=95.25, \
                    llcrnrlat=37, \
                    urcrnrlon=97.25, \
                    urcrnrlat=38, \
                    resolution='h',epsg = 4269)
            parallels = np.arange(37,38,1.)
            meridians = np.arange(95.5,97.5,0.5)
            xpixels = 1000
           # map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)
            eastings, northings = map(eastings, northings)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[:,1], data[:,0])
            mins = np.max(data[:,2])
            #data_int[data_int<np.max(data_int)*0.000001]=np.nan
            import matplotlib.colors as colors
            import matplotlib.tri as tri
            #mask = np.ma.masked_where(data_int < 0.4, data_int)
            #mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang = tri.Triangulation(x, y)
            isbad = np.less(data_int, 0.01)
            mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
            triang.set_mask(mask)
            plt.tricontourf(triang, data_int, cmap='YlOrRd')
            #plt.tricontourf(x,y, data_int, cmap='hot',alpha=0.6)

            #plt.tricontourf(x,y, data_int, cmap='hot',norm=colors.Normalize(vmin=0.1, vmax=1.1))
            plt.colorbar()
            plt.title(path_in_str)
            ax = plt.gca()
            np1 = [116, 61, 91]
            x, y = map(96.476,37.529)

            beach1 = beach(np1, xy=(x, y), width=0.05)
            ax.add_collection(beach1)
            xpixels = 1000
            map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

            plt.show()


if len(sys.argv)<3:
    print("input: eventname plot_name,\
     available plot_name: movie, sembmax, semblance, interactive_max, cluster")
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
        plot_scatter()
