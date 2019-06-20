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
from matplotlib.patches import Circle, Polygon
w=25480390.0

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



stats = 'event.stations'

f = open("event.stations") # open your file

f=open('event.stations',"r")
lines=f.readlines()
result=[]
for x in lines:
    result.append(x.split(' ')[:])
f.close()
event_cor_y = 37.6199
event_cor_x = 95.8853
lat_0 = 37.6199
lon_0 = 95.8853

map = Basemap(width=29000000,height=29000000,
            resolution='l',projection='stere',\
            lat_ts=event_cor_y,lat_0=event_cor_y,lon_0=event_cor_x)
map.drawcoastlines()

ax = plt.gca()
np1 = [101, 60, 83]
x, y = map(95.8853,37.6199)

beach1 = beach(np1, xy=(x, y), width=700000)
ax.add_collection(beach1)

x = []
y = []
z = []
for i in result:
    x.append(i[:][2])
    y.append(i[:][1])
    z.append(int(i[:][3]))


colors = cm.nipy_spectral(np.linspace(0,1,np.max(z)+1))

x = np.asarray(x)
y = np.asarray(y)

x, y = map(x,y)
map.scatter(x,y,c=[colors[index] for index in z])
try:
    plt.text(x[0],y[0],'r'+str(data[0,0])[:], fontsize=12)
except:
    pass

radii = [1000,5000,10000]
for radius in radii:
    equi(map, lon_0, lat_0, radius,lw=2.)
plt.show()
