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

def plot_movie():
    if len(sys.argv)<4:
        print "missing input arrayname"
    else:
        if sys.argv[3] == 'combined':
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
        else:
            rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/' + str(sys.argv[3])
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
            plt.savefig(path_in_str+'.pdf', bbox_inches='tight')

def plot_sembmax():
    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    data = num.loadtxt(rel+'sembmax_0.txt', delimiter=' ', skiprows=5)
    eastings = data[:,2]
    northings =  data[:,1]

    map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
            resolution='h',epsg = 4269)

    X,Y = np.meshgrid(eastings, northings)

    eastings, northings = map(X, Y)
    map.drawcoastlines(color='b',linewidth=1)

    x, y = map(data[:,2], data[:,1])
    l = range(0,num.shape(data[:,2])[0])

    ps = map.scatter(x,y,marker='o',c=l, s=data[:,3]*8000, cmap='seismic', vmin= num.max(data[:,3])*0.66)
    xpixels = 1000
    map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

    cbar = map.colorbar(ps,location='bottom',pad="5%", label='Time [s]')
    plt.savefig(rel+'semblance_max_0.pdf', bbox_inches='tight')
    plt.show()
    try:
        rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
        data = num.loadtxt(rel+'sembmax_1.txt', delimiter=' ', skiprows=5)
        eastings = data[:,2]
        northings =  data[:,1]

        map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
                resolution='h',epsg = 4269)

        X,Y = np.meshgrid(eastings, northings)

        eastings, northings = map(X, Y)
        map.drawcoastlines(color='b',linewidth=1)

        x, y = map(data[:,2], data[:,1])
        l = range(0,num.shape(data[:,2])[0])

        ps = map.scatter(x,y,marker='o',c=l, s=data[:,3]*8000, cmap='seismic', vmin= num.max(data[:,3])*0.66)
        xpixels = 1000
        map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

        cbar = map.colorbar(ps,location='bottom',pad="5%", label='Time [s]')
        plt.savefig(rel+'semblance_max_1.pdf', bbox_inches='tight')
        plt.show()
    except:
        pass


def plot_movingsembmax():
    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    data = num.loadtxt(rel+'sembmax_0.txt', delimiter=' ', skiprows=5)
    eastings = data[:,2]
    northings =  data[:,1]
    xpixels = 1000
    map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
            resolution='h',epsg = 4269)

    X,Y = np.meshgrid(eastings, northings)

    eastings, northings = map(X, Y)
    map.drawcoastlines(color='b',linewidth=1)
    map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

    x, y = map(data[:,2], data[:,1])
    size = num.shape(data[:,2])[0]
    l = range(0,size)

    scat = map.scatter(x,y,marker='o',c=l, cmap='jet', s=data[:,3]*10900)
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
        data = num.loadtxt(rel+'sembmax_1.txt', delimiter=' ', skiprows=5)
        eastings = data[:,2]
        northings =  data[:,1]
        xpixels = 1000
        map = Basemap(projection='merc', llcrnrlon=num.min(eastings),llcrnrlat=num.min(northings),urcrnrlon=num.max(eastings),urcrnrlat=num.max(northings),
                resolution='h',epsg = 4269)

        X,Y = np.meshgrid(eastings, northings)

        eastings, northings = map(X, Y)
        map.drawcoastlines(color='b',linewidth=1)
        map.arcgisimage(service='World_Shaded_Relief', xpixels = xpixels, verbose= False)

        x, y = map(data[:,2], data[:,1])
        size = num.shape(data[:,2])[0]
        l = range(0,size)

        scat = map.scatter(x,y,marker='o',c=l, cmap='jet', s=data[:,3]*10900)
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

    rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
    astf = num.loadtxt(rel+'sembmax_0.txt', delimiter=' ', skiprows=5)
    astf_data= astf[:,3]

    fig = plt.figure(figsize=(15,15))

    plt.plot(astf_data)
    plt.ylabel('Beampower')

    plt.xlabel('Time [s]')

    plt.savefig(rel+'semblance_0.pdf', bbox_inches='tight')
    plt.show()
    try:
        rel = 'events/'+ str(sys.argv[1]) + '/work/semblance/'
        astf = num.loadtxt(rel+'sembmax_1.txt', delimiter=' ', skiprows=5)
        fig = plt.figure(figsize=(15,15))

        plt.plot(astf_data)
        plt.ylabel('Beampower')

        plt.xlabel('Time [s]')

        plt.savefig(rel+'semblance_1.pdf', bbox_inches='tight')
        plt.show()
    except:
        pass

if len(sys.argv)<3:
    print "input: eventname plot_name"

event = sys.argv[1]
if sys.argv[2] == 'plot_movie':
    plot_movie()
elif sys.argv[2] == 'plot_sembmax':
    plot_sembmax()
elif sys.argv[2] == 'plot_semb':
    plot_semb()
elif sys.argv[2] == 'plot_movingmax':
    plot_movingsembmax()
