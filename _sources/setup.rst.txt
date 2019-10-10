
Setup
========


Installation from source
^^^^^^^^^^^^^^^^^^^^^^^^

'git clone https://git.pyrocko.org/Palantiri'
and than in the cloned directory::

  sudo python3 setup.py install


Note: This tool is written in python3.
A python2 branch is available as backup but not maintained.


Prerequisites
^^^^^^^^^^^^^

It is necessary to install the pyrocko <https://github.com/pyrocko/pyrocko>`_. software package, which is used to handle
the basic waveform handling and manipulation. Please follow the pyrocko packages installation guide, before installing this software.

All prerequisites listed for the pyrocko software are assumed to be installed for the usage of this software package.

For some advanced functionality (e.g. array response analysis) it is also necessary to install the obspy package:
obspy: <https://github.com/obspy/obspy>`_
Some further basic requirements for plotting can be installed with::

  sudo pip3 install pyproj basemap affine

These are only necessary for plotting and can be omitted.
