
Processing
==========

After step 2 you will have a folder named after the specific event in the events subfolder. This
eventfolder will contain all work and data specific to this event.




processing steps:

step 1)

		palantiri_eventsearch: - search for earthquakes to process in global catalogs
                                    - searchparameter are defined in global.conf
                                    - possible parameter: - date from to
                                                     - magnitude
                                                     - catalog
                                                     - number of result


step 2)

		palantiri_create eventid: - creates eventfolder in events
                                            - use eventid from search to start creating event directory structure
                                                (copies example.config from skeleton folder)


For step 3) three options exist:

a) Download real data with the download tool (faster, less stations):

		palantiri_down

b) use synthetics but with a station distributions from a real data case:

		For this you will need a greensfunctions store that is pre-calculated with the fomosto tool from pyrocko (https://pyrocko.org/docs/current/apps/fomosto/index.html).
		Several already pre-calculated stores for download can be found at http://kinherd.org:8080/gfws/static/stores/
		This possibilty assumes also that you downloaded data with a) or b), as the real station distributions will be used for the synthetics.
		Please make sure to set the pyrocko_download option in the event config file to 1 if you downloaded data with this command.
		Also the noise of the real data will be used to pertub the synthetics if you select the option in the event config file.


c) use synthetics from a scenario generator:

		You can also use the output of the pyrocko scenario generator colosseo.
		After you followed the steps to create a purely synthetic dataset at https://pyrocko.org/docs/current/apps/colosseo/index.html
		you have to give the scenario.yml of the scenario as input in the event configuration file and set the variable colosseo input to
		1. Please make sure that you unset other types of input. Also give the greensfunctions store used in the synthetic input file
		(eventname.syn). Disregard all other parameters in the synthetic input file for this case, as the generated event from the scenario
		will be used. You will need to merge all mseed files in the scenarios waveform directory into one file called scenario.mseed, located
		in the same folder as the scenario.yml. This can be done with jackseis or simply by using cat.

step 4) Clustering the stations into virtual arrays, based on the input in the eventname.config. This is handled as an optimization problem and returns a set of virtual arrays.

		palantiri_cluster eventname: - clustering of stations to automatically build virtual arrays
                                                     - configuration parameter in the event config file (eventname.config)
                                                     - possible parameter:


STEP 5) The last step is the actual processing.

The next steps are based on the input you have chosen before. Be sure to not mix different types of input. Remove or move the folders eventname/cluster and
eventname/work if you want to start over for different input or array setup.
Again be careful to check the eventname.config file in the event folder and adjust it your liking.
Note that the option for several filters are build in. With this option the processing will be done separately for different filter setups
and according outputs are generated. An arbitrary number of filter can used. The filter parameters names are assumed to be consecutively numbered.   This processing of different filter settings is useful for exploring e.g. high- and low-frequency content.
Also note that several depths can be selected to iterate over. Else only one planar equi-depth grid is considered for the semblance and traveltime
calculation. If several depths are chosen the processing will be repeated for each depth and the semblance will be output for each depth.
Arrays can be weighted by pre-event noise variance and azimuthal coverage.


The semblance output is located in the eventfolder/work/semblance as txt files with the ending ".asc". They are

	First the data of each array can be cross-correlated. Stations under the threshold (xcorrtreshold) given in the eventname.config are disregarded. They are crosscorrelated by default to the first station of the array but a reference station can be manually given to each array. 	xcorr=1 enables a correction of timeshifts at each based on cross correlations. If also autoxcorrcorrectur = 1 is selected for each array a manual picking of phase onsets is done before the processing. This will return a reference waveform of one of the stations
	in the virtual array in a figure and a snuffler window.  Marker for STA/LTA and theoretical phase onsets will be given.
	After closing both figures, the user can then input a manual traveltime shift in second in regard to the xcorr window start (also markers in the snuffler). The traveltimes for this array will than be statically corrected using this manual selected value. Both methods allows for handling of velocity model inadequacies.


	Second the traveltimes for each gridpoint to each station will be pre-calculated. This can take some time, depending on your gridsize. The traveltime grids are
	saved automatically in the folder tttgrids for each array separately. They will automatically be loaded in when starting step 5 again. This is very useful for synthetic
	test as it saves a lot of time. If you change the setup of arrays however you will have to delete the saved tttgrid files for the affected arrays. If the dimensions of the grid change they will have to be calculated again as well.

	Lastly for the semblance calculation two options exists. Firstly the semblance can be calculated for each array separately and then combined. The combination can be weighted by the average SNR of the arrays if the option
	is choosen in the eventname.config. The output are grids for each timestep of semblance which are stored in eventname/work/semblance for each array in a different folder with the ending
	.asc. The combined semblance for all arrays can be found directly in eventname/work/semblance also with the ending .asc. If you used multiple filter, the files will have a numeral matching the
	listing of the filter. Also for each depth chosen a different output will be generated.

The second option will combine all stations from all arrays into


		bat process eventname: - ARRAYPROCESSING OF THE EVENT
                                                     - CONFIGURATION PARAMETER IN THE EVENT CONFIG FILE
