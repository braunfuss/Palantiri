
Palantiri Manual
================

Palantiri is an console command based open source toolbox for seismological backprojection to investigate source properties
based on teleseismic data. Palantiri allows for downloading data and clustering stations into synthetics arrays.
Bootstrapping of the weight of this arrays allow to investigate the uncertainty of the backprojection results.
The tool allows for a number of fast synthetic tests, which are made possible using Greensfunctions stores from the Pyrocko package.
For fast processing traveltime grids are pre-calculated. Arbitrary number of processes for different filter settings can be run.
Backprojections can be further carried out on a single planar grid or in 3-D.



.. raw:: html

   <div style="clear:both"></div>

Features
--------

.. raw:: html

   <div class="keywords">
   <span>time domain backprojection</span>
   <span>array forming by clustering</span>
    <span>automatic download of data </span>
   </div>

Contents
--------

.. toctree::
   :maxdepth: 2

   overview
   setup
   configuration/index
   processing

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Literature
------------------
  KRÜGER, Frank; OHRNBERGER, Matthias.
  Tracking the rupture of the M w= 9.3 Sumatra earthquake over 1,150 km at teleseismic distance.
  Nature, 2005, 435. Jg., Nr. 7044, S. 937.

  Rössler, D., Krueger, F., Ohrnberger, M., & Ehlert, L. (2010).
  Rapid characterisation of large earthquakes by multiple seismic broadband arrays.
  Natural Hazards and Earth System Sciences, 10(4), 923-932.


More resources:
---------------

  * `Green's Mill - Online resource for pre-calculated Green's functions <https://greens-mill.pyrocko.org>`_
