# seismoBat
### _A seismological backprojection tool


## Installation



### Installation from source

```
git clone https://github.org/braunfuss/seismoBat

Prequsites:
pyrocko: https://github.com/pyrocko/pyrocko
obspy: https://github.com/obspy/obspy

Note: This tool is still writen in python2, therefore you
need to install obspy and pyrocko with python2!
Please make sure that current version of both obspy and pyrocko
are installed. 
```


## Quick Manual


- for general options, regarding data and event selection choose parameters in global.conf, which is in the seismoBat main folder.
- event dependent options you can change in the config file of the event in the eventfolder, which will be created after step 2



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
		
	
		python arraytool.py pyrocko_download eventfoldername

b) Download  real data with obspy (slower, more stations). For this three commands are needed, first:
		

		1) python arraytool.py getstations eventfoldername: - search for station in iris and webdc with searchparameters defined in global.conf
		                                                 - possible parameter: - mindistance from event
		                                                                    - maxdistance from event
		                                                                    - networks for blacklist

  		2) python arraytool.py getdata eventfoldername: - acquisition of archived waveformdata from iris and eida saved as sds in the eventfolder
		                                             - copy keyfoldername from step 3 into global.conf
		                                             - possible parameter: - keyfolder
		                                             - duration of waveformdata


	       3) python arraytool.py getmeta eventfoldername: - create file with metainformation from every station which has data in sds structure


c) use synthetics but station distributions from a real data case: 
		For this you will need a greensfunctions store that is pre-calculated with the fomosto tool from pyrocko (https://pyrocko.org/docs/current/apps/fomosto/index.html).
		Several already pre-calculated stores for download can be found at http://kinherd.org:8080/gfws/static/stores/
		This possibilty assumes also that you downloaded data with a) or b).


d) use synthetics from a scenario generator:

		You can also use the output of the pyrocko scenario generator colosseo.
		After you followed the steps to create a purely synthetic dataset at https://pyrocko.org/docs/current/apps/colosseo/index.html
		you can give the scenario.yml to 
		

step 4) 

		python arraytool.py cluster eventfoldername: - clustering of stations to automatically build synthetic arrays (gives you arrayconfiguration + overview plot)
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

STEP 5) 

		python arraytool.py process eventfoldername: - ARRAYPROCESSING OF THE EVENT
                                                     - CONFIGURATION PARAMETER IN THE EVENT CONFIG FILE
						     - this can also be done with synthetic input data

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
 24118 Kiel, Germanydam, Germany

```

