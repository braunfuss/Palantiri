
Setup
========


Installation from source
^^^^^^^^^^^^^^^^^^^^^^^^


git clone https://git.pyrocko.org/Palantiri
and than in the cloned directory:
sudo python3 setup.py install


Note: This tool is written python3.
A python2 branch is available as backup but not maintained.


Prerequisites
^^^^^^^^^^^^^

It is necessary to install the pyrocko software package, which is used to handle
the basic waveform handling and manipulation. Please see that packages installation guide on:
pyrocko: https://git.pyrocko.org/pyrocko/
All prerequisites listed for the pyrocko software are assumed to be installed for the usage of this software package.

For some advanced functionality (e.g. array response analysis) it is also necessary to install the obspy package:
obspy: https://github.com/obspy/obspy

Some further basic installations which are needed are:
sudo pip3 install pyproj basemap affine

basemap is only necessary for plotting and can be omitted.
