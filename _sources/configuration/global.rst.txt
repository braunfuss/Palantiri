global.conf
==============

First the file global.conf, which is located in the palantiri main workdir gives options for general download options, regarding data and event selection.
An example file will be created when using the palantiri_init command. The user can give a date_min, date_max and a minimum magnitude to specify for which events should be searched.
All FDSN catalogs are supported.

An example configuration file::

  [eventsearchparameter]

  date_min = 1995-11-22T20:00:00.0
  date_max = 1995-11-22T24:00:00.0
  magmin = 2.5
  catalog =  GCMT
  resultlimit = 10

  [stationsearchparameter]

  blacklist = SY

In this example we select [eventsearchparameter] for an earthquake above magnitude 2.5, occurring between 1995-11-22T20:00:00.0 and 1995-11-22T24:00:00.0 in the GCMT catalog.
We display only the first 10 results.

The stationsearchparameter 'blacklist' is used in the clustering to ignore stations in all processing.
A comma separated list is expected.
