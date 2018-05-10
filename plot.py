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



def plot_movie():
    if len(sys.argv)<3:
        print "missing input arrayname"
    else:
        rel = 'events'+ str(sys.argv[0]) + 'work/semblance/' + str(sys.argv[3])
        pathlist = Path(rel).glob('**/*.ASC')
        for path in sorted(pathlist):
            path_in_str = str(path)
            data = num.loadtxt(path_in_str, delimiter=' ', skiprows=5)

            eastings = data[:,1]
            northings =  data[:,0]
            plt.figure()
            map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
                    resolution='h')
            parallels = np.arange(num.min(northings),num.max(northings),1.)
            meridians = np.arange(num.min(eastings),num.max(eastings),1.)

            eastings, northings = map(eastings, northings)
            map.drawcoastlines(color='b',linewidth=3)
            map.drawparallels(parallels,labels=[1,0,0,0],fontsize=22)
            map.drawmeridians(meridians,labels=[1,1,0,1],fontsize=22)
            x, y = map(data[::5,1], data[::5,0])
            mins = np.max(data[:,3])
            plt.tricontourf(x,y, data[::5,3], vmin=mins*0.6)
            plt.title(path_in_str)
            plt.savefig(rel+path_in_str+'.pdf', bbox_inches='tight')

def plot_sembmax():
    rel = 'events'+ str(sys.argv[0]) + 'work/semblance/sembmax_1.txt'
    astf = num.loadtxt(rel, delimiter=' ', skiprows=5)
    data= astf[:,3]
    eastings = data[:,2]
    northings =  data[:,1]

    map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
            resolution='h',epsg = 4269)

    X,Y = np.meshgrid(eastings, northings)

    eastings, northings = map(X, Y)
    map.drawcoastlines(color='b',linewidth=1)

    x, y = map(data[0:time,2], data[0:time,1])
    #xs, ys = map(74.219,39.226)
    l = range(0,num.shape(data[0:time,2])[0])

    ps = map.scatter(x,y,marker='o',c=l, s=data[:,3]*8000, cmap='seismic', vmin= num.max(data[:,3])*0.66)

    #map.scatter(xs,ys, marker='o', linewidth=6)
    xpixels = 1000
    map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

    cbar = map.colorbar(ps,location='bottom',pad="5%", label='Time [s]')
    plt.savefig(rel+'semblance_max_1.pdf', bbox_inches='tight')
    try:
        rel = 'events'+ str(sys.argv[0]) + 'work/semblance/sembmax_2.txt'
        astf = num.loadtxt(rel, delimiter=' ', skiprows=5)
        data= astf[:,3]
        eastings = data[:,2]
        northings =  data[:,1]

        map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
                resolution='h',epsg = 4269)

        X,Y = np.meshgrid(eastings, northings)

        eastings, northings = map(X, Y)
        map.drawcoastlines(color='b',linewidth=1)

        x, y = map(data[0:time,2], data[0:time,1])
        #xs, ys = map(74.219,39.226)
        l = range(0,num.shape(data[0:time,2])[0])

        ps = map.scatter(x,y,marker='o',c=l, s=data[:,3]*8000, cmap='seismic', vmin= num.max(data[:,3])*0.66)

        #map.scatter(xs,ys, marker='o', linewidth=6)
        xpixels = 1000
        map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

        cbar = map.colorbar(ps,location='bottom',pad="5%", label='Time [s]')
        plt.savefig(rel+'semblance_max_2.pdf', bbox_inches='tight')
    except:
        pass


def plot_semb():

    rel = 'events'+ str(sys.argv[0]) + 'work/semblance/sembmax_1.txt'
    astf = num.loadtxt(rel, delimiter=' ', skiprows=5)
    astf_data= astf[:,3]

    fig = plt.figure(figsize=(15,15))

    plt.plot(astf_data)
    plt.ylabel('Beampower')

    plt.xlabel('Time [s]')

    plt.savefig(rel+'semblance_1.pdf', bbox_inches='tight')
    try:
        rel = 'events'+ str(sys.argv[0]) + 'work/semblance/sembmax_2.txt'
        astf = num.loadtxt(rel, delimiter=' ', skiprows=5)
        fig = plt.figure(figsize=(15,15))

        plt.plot(astf_data)
        plt.ylabel('Beampower')

        plt.xlabel('Time [s]')

        plt.savefig(rel+'semblance_2.pdf', bbox_inches='tight')
    except:
        pass

if len(sys.argv)<2:
    print "input: eventname plot_name"

event = sys.argv[0]

if sys.argv[1] == 'plot_movie'
    plot_movie()
elif: sys.argv[1] == 'plot_sembmax'
    plot_sembmax()
elif: sys.argv[1] == 'plot_semb'
    plot_semb()
