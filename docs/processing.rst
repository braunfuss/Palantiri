Pre-processing
==============

Please make sure that you have configured the configuration files beforehand, as detailed in the configuration section.
Palantiri is a command line tool.

As a first step we need to search for an event to work on:

1. Event search
---------------

The event search is carried out by:

.. highlight:: console

::

    $ palantiri_eventsearch
    Usage: palantiri_eventsearch


A list is returned in the terminal with the ISC event ids of the found earthquakes.

2. Event folder creation
------------------------

 We can now use the id found in step 1 to create an event folder for the selected earthquake by:

.. highlight:: console

::

    $ palantiri_create
    Usage: palantiri_create eventid


An eventfolder is created in the palantiri main workdir under the events folder with the event name and date as folder name. The folder contains some example configuration files
and a basic folder structure::

    eventname_date/
    |-- eventname_date.config 		# (1)
		|-- eventname_date.syn			  # (2)
		|--eventname_date.origin   		# (3)
    `-- data/
		`-- work/
			`-- semblance


2 a) Data download
---------------------

Download real data with the download tool:

.. highlight:: console

::

		$ palantiri_down
    usage: palantiri_down [options] <eventname_date>  "<YYYY-MM-DD HH:MM:SS>" <lat> <lon> \\
                                   <depth_km> <radius_km> <fmin_hz> \\
                                   <sampling_rate_hz> \\
                                   <eventname> [--]

           palantiri_down [options] <eventname_date> "<YYYY-MM-DD HH:MM:SS>" <radius_km> <fmin_hz> \\
                                   <sampling_rate_hz> <eventname> [--]

           palantiri_down [options] <eventname_date>  <catalog-eventname> <radius_km> <fmin_hz> \\
                                   <sampling_rate_hz> <eventname> [--]

           palantiri_down [options] --window="<YYYY-MM-DD HH:MM:SS, YYYY-MM-DD HH:MM:SS>" \\
                                   <eventname_date>  <lat> <lon> <radius_km> <fmin_hz> \\
                                   <sampling_rate_hz> <eventname> [--]

2 b) Synthetic data generation
-------------------------------

For generating synthetic data for palantiri, two options exist.
First you can use the station distributions you obtained from downloading using the palantiri_down command.
This is recommend to check the array response and the sensitivity in preparation for a real data analysis.

For this you will need a greensfunctions store that is pre-calculated with the fomosto tool from pyrocko (https://pyrocko.org/docs/current/apps/fomosto/index.html).
Several already pre-calculated stores for download can be found at http://kinherd.org:8080/gfws/static/stores/
This possibilty assumes also that you downloaded data with a) or b), as the real station distributions will be used for the synthetics.
Please make sure to set the pyrocko_download option in the event config file to 1 if you downloaded data with this command.
Also the noise of the real data will be used to perturb the synthetics if you select the option in the event config file.

As a second option we offer to use the output of the colosseo scneario generator.

You can also use the output of the pyrocko scenario generator colosseo.
After you followed the steps to create a purely synthetic dataset at https://pyrocko.org/docs/current/apps/colosseo/index.html
you have to give the scenario.yml of the scenario as input in the event configuration file and set the variable colosseo input to
1. Please make sure that you unset other types of input. Also give the greensfunctions store used in the synthetic input file
(eventname_date.syn). Disregard all other parameters in the synthetic input file for this case, as the generated event from the scenario
will be used. You will need to merge all mseed files in the scenarios waveform directory into one file called scenario.mseed, located
in the same folder as the scenario.yml. This can be done with jackseis or simply by using cat.


3. Array clustering
-------------------------------

In this step we cluster the stations into virtual arrays, based on the configuration in the eventname_date.config. This is handled as an kmeans problem and returns a set of virtual arrays.
The command to cluster the stations into virtual arrays is:

.. highlight:: console

::

		$ palantiri_cluster
		Usage: palantiri_cluster eventname_date



The desired stations in each array can be given/modified in the eventname_date.config file, also allowing for manual array creation.
As input a comma separated list of station names is expected in the format::

  Network.Station..Channel

The output arrays are numbered and assigned a prefix ``r``. Note that the numbering might not be consecutive as some arrays will be disregarded after the clustering because of the settings in the configuration file (distance from source to stations, distance between arrays).

Each array can also be given a reference station (``arrayrefstation``) which will be used for cross-correlation and plotting purposes.
If not given the station which is closest to the center of the array will be used as reference.


Processing
==========
The last step is the actual processing.
This chapter describes the main processing. After the pre-processing you will have a folder named after the specific event in the events subfolder and your eventname_date.config file contains a list of arrays.
The eventfolder will contain all work and data specific to this event. If you reprocess a certain event the results will be overwritten.
For beamforming several methods are incorporated, including linear, phase-weighted and coherence based stacking.

The MUSIC algorithm is at this stage partly supported but will be fully implemented in later versions.


The next steps are based on the input you have chosen before. Be sure to not mix different types of input. Remove or move the folders eventname_date/cluster and
eventname_date/work if you want to start over for different input or array setup.
Again be careful to check the eventname_date.config file in the event folder and adjust it your liking.
Note that the option for several filters are build in. With this option the processing will be done separately for different filter setups
and according outputs are generated. An arbitrary number of filter can used. The filter parameters names are assumed to be consecutively numbered.   This processing of different filter settings is useful for exploring e.g. high- and low-frequency content.
Also note that several depths can be selected to iterate over. Else only one -planar (equi-depth)- grid is considered for the semblance and traveltime
calculation. If several depths are chosen the processing will be repeated for each depth and the semblance will be output for each depth.
Arrays can be weighted by pre-event noise variance and azimuthal coverage.


The semblance output is located in the eventfolder/work/semblance as txt files with the ending ``.asc``. They are

First the data of each array can be cross-correlated. Stations under the threshold (xcorrtreshold) given in the eventname_date.config are disregarded. They are crosscorrelated by default to the first station of the array but a reference station can be manually given to each array. 	xcorr=1 enables a correction of timeshifts at each based on cross correlations. If also autoxcorrcorrectur = 1 is selected for each array a manual picking of phase onsets is done before the processing. This will return a reference waveform of one of the stations
in the virtual array in a figure and a snuffler window.  Marker for STA/LTA and theoretical phase onsets will be given.
After closing both figures, the user can then input a manual traveltime shift in second in regard to the xcorr window start (also markers in the snuffler). The traveltimes for this array will than be statically corrected using this manual selected value. Both methods allows for handling of velocity model inadequacies.


Second the traveltimes for each gridpoint to each station will be pre-calculated. This can take some time, depending on your gridsize. The traveltime grids are
saved automatically in the folder tttgrids for each array separately. They will automatically be loaded in when starting step 5 again. This is very useful for synthetic
test as it saves a lot of time. If you change the setup of arrays however you will have to delete the saved tttgrid files for the affected arrays. If the dimensions of the grid change they will have to be calculated again as well.

Lastly for the semblance calculation two options exists. Firstly the semblance can be calculated for each array separately and then combined. The combination can be weighted by the average SNR of the arrays if the option
is chosen in the eventname_date.config. The output are grids for each timestep of semblance which are stored in eventname_date/work/semblance for each array in a different folder with the ending
.asc. The combined semblance for all arrays can be found directly in eventname_date/work/semblance also with the ending ``*.asc``. If you used multiple filter, the files will have a numeral matching the
listing of the filter. Also for each depth chosen a different output will be generated.

The actual processing is carried out by calling the bat command:

.. highlight:: console

::

		$ bat
		Usage: bat process eventname_date
