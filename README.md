# seismoBat
### A seismological backprojection array tool


## Installation



### Installation from source

```
git clone https://github.org/braunfuss/seismoBat

Prequsites:
pyrocko: https://github.com/pyrocko/pyrocko
obspy: https://github.com/obspy/obspy

Basics that need to be installed:

sudo pip install pyproj scipy matplotlib numpy basemap affine 

basemap is only necessary for plotting and can be omitted

Please make sure that current version of both obspy and pyrocko
are installed. 

Note: This tool is writen in python2, therefore you
need to install obspy and pyrocko with python2!

```


## Quick Manual

This tool is intended to 
Before processing, this tool can download waveform data and cluster them automatically into virtual arrays.
Also purely synthetic tests are possible, either with synthetics generated for stations of a real case or
with input from the pyrocko colosseo scenario generator (https://pyrocko.org/docs/current/apps/colosseo/index.html).

There are three config files which specify the user input:

First, for general options, regarding data and event selection choose parameters in global.conf, which is in the seismoBat main folder.
This options are used for steps 0-2. After step 2 you will have a folder named after the specifc event in the events subfolder. This
eventfolder will contain all work and data specifc to this event. 

Secondly, event dependent options can be changed in the config file (eventname.config) of the event in the eventfolder, which will be created after step 2. 
Please make sure to investigate this configuration file closely before further processing. 

And lastly, synthetic sources can be specified in the syn config file (eventname.syn) in the eventfolder, also created after step 2. A user definded number of
RectangularSources or DoubleCouble sources can used to create synthetic data with free user input of the source parameters. 
See the pyrocko documentation for details about the sources. This configuration file also contains paths to greensfunction stores.



processing steps:

step 0) 

		python arraytool.py list: - lists all events in the eventfolder (config and orig file must exist) (optional, just to see which events are available to process)

step 1) 

		python arraytool.py search: - search for earthquakes to process in global catalogs
                                    - searchparameter are defined in global.conf
                                    - possible parameter: - date from to
                                                     - magnitude
                                                     - catalog
                                                     - number of result

                            
step 2) 

		python arraytool.py create eventid: - creates eventfolder in events 
                                            - use eventid from search to start creating event directory structure
                                                (copies example.config from skeleton folder)
                                    

For step 3) four options exsist:
	
a) Download real data with pyrocko (faster, less stations):
		
	
		python arraytool.py pyrocko_download eventname

b) Download  real data with obspy (slower, more stations). For this three commands are needed, first:
		

		1) python arraytool.py getstations eventname: - search for station in iris and webdc with searchparameters defined in global.conf
		                                                 - possible parameter: - mindistance from event
		                                                                    - maxdistance from event
		                                                                    - networks for blacklist

		2) python arraytool.py getdata eventname: - acquisition of archived waveformdata from iris and eida saved as sds in the eventfolder
		                                             - copy keyfoldername from step 3 into global.conf
		                                             - possible parameter: - keyfolder
		                                             - duration of waveformdata


		3) python arraytool.py getmeta eventname: - create file with metainformation from every station which has data in sds structure


c) use synthetics but station distributions from a real data case: 

		For this you will need a greensfunctions store that is pre-calculated with the fomosto tool from pyrocko (https://pyrocko.org/docs/current/apps/fomosto/index.html).
		Several already pre-calculated stores for download can be found at http://kinherd.org:8080/gfws/static/stores/
		This possibilty assumes also that you downloaded data with a) or b), as the real station distributions will be used for the synthetics.
		Please make sure to set the pyrocko_download option in the event config file to 1 if you downloaded data with this command.
		Also the noise of the real data will be used to pertub the synthetics if you select the option in the event config file.


d) use synthetics from a scenario generator:

		You can also use the output of the pyrocko scenario generator colosseo.
		After you followed the steps to create a purely synthetic dataset at https://pyrocko.org/docs/current/apps/colosseo/index.html
		you have to give the scenario.yml of the scenario as input in the event configuration file and set the variable colosseo input to
		1. Please make sure that you unset other types of input. Also give the greensfunctions store used in the synthetic input file 
		(eventname.syn). Disregard all other parameters in the synthetic input file for this case, as the generated event from the scenario
		will be used. 
		

The next steps are based on the input you have choosen before. Be sure to not mix different types of input. Remove or move the folders eventname/cluster and
eventname/work if you want to start over for different input or array setup. 
Again be careful to check the eventname.config file and adjust it your liking.
Note that the option for several filters are build in. With this option the processing will be done seperatly for different filter setups
and according outputs are generated. This is useful for divding into high- and low-frequency content. 
Also note that several depths can be selected to iterate over. Right now for one run only planar equi-depth grid are considered for the semblance
calculation. If several depths are choosen the processing will be repeted for each depth.
step 4) Clustering the stations into virtual arrays, based on the input in the eventname.config. This is handled as an optimization problem
	and returns the best solution.

		python arraytool.py cluster eventname: - clustering of stations to automatically build virtual arrays 
                                                     - configuration parameter in the event config file (eventname.config)
                                                     - possible parameter:  - clustermethod (distance zoning/ kmeans)
                                                                              distance zoning:
                                                                                    beginminstationdistance = 1
                                                                                    endminstationdistance=59.9
                                                                                    beginmaxstationdistance=60
                                                                                    endmaxstationdistance=90
                                                            
                                                      kmeans:
                                                           maxcluster = 2 (amount of clusters you would like to have at the end of clustering)
                                                           minstationaroundinitialcluster = 5 (minimum station around initial cluster center to be one inital cluster)
                                                           initialstationdistance = 10 (stations must be 10 degree around initial cluster center)
                                                           cutoff = 30              (if no final result for kmean then run only 30 times and take last result)
                                                           runs = 5                 (repeatings for kmean clustering to get best result)
                                                           centroidmindistance = 20 (minimum distance of initial centroids,unit in degree)
                                                           comparedelta = 2       
                                                           stationdistance = 10   (maximum distance from station to cluster center)
                                                           minclusterstation = 10 (minimum stations per cluster)

STEP 5) The last step is the actual processing. 
	First the data of each array will be cross-correlated. Stations under the threshold (xcorrtreshold) given in the eventname.config
	are disregarded. If autoxcorrcorrectur = 1 is selected for each array a manual qualtiy control is done before and will return a reference waveform of one of the stations
	in the virtual array. A figure and a snuffler window will open for the user to investigate the traveltime. Marker for STA/LTA and theortical phase onsets will be given.
	After closing both figures, the user can then input a manual traveltime shift in second in regard to the xcorr window start (also markers in the snuffler).
	

	Second the traveltimes for each gridpoint to each station will be calculated. This can take some time, depending on your gridsize. Therefor the traveltime grids are 
	saved automatically in the folder tttgrids for each array seperatly. They will automatically be loaded in when starting step 5 again. This is very useful for synthetic
	test as it saves a lot of time. If you change the setup of the arrays however (with step 4) you will have to delete the saved tttgrid files for the affected arrays.

	Lastly the semblance will be calculated first for each array seperatly and then combined. The combination can be weighted by the average SNR of the arrays if the option
	is choosen in the eventname.config. The output are grids for each timestep of semblance which are stored in eventname/work/semblance for each array in a different folder with the ending
	.asc. The combined semblance for all arrays can be found directly in eventname/work/semblance also with the ending .asc. If you used multiple filter, the files will have a numeral matching the 
	listing of the filter. Also for each depth choosen a different output will be generated.
	


		python arraytool.py process eventname: - ARRAYPROCESSING OF THE EVENT
                                                     - CONFIGURATION PARAMETER IN THE EVENT CONFIG FILE
						     

## Documentation

Documentation and usage examples are available online soon


## Citation


## License 
GNU General Public License, Version 3, 29 June 2007

Copyright © 2018 University Potsdam, Potsdam, Germany and  University of Kiel, Kiel, Germany

seismoBat is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
seismoBat is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

## Contact
* Andreas Steinberg; 
  andreas.steinberg@ifg.uni-kiel.de

* Frank Krüger; 
  kruegerf@geo.uni-potsdam.de


```
 University of Kiel
 Institute of Geosciences
 Otto-Hahn-Platz 1
 24118 Kiel, Germany, Germany

```

