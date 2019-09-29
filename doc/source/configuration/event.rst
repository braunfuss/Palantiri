event.conf
===========

Secondly, event dependent options can be changed in the config file (eventname.config) of the event in the eventfolder, which will be created after step 2.
Please make sure to investigate this configuration file closely before further processing.

And lastly, synthetic sources can be specified in the syn config file (eventname.syn) in the eventfolder, also created after step 2. A user definded number of
RectangularSources or DoubleCouble sources can used to create synthetic data with free user input of the source parameters.
See the pyrocko documentation for details about the sources. This configuration file also contains paths to greensfunction stores.

The main event configuration file is created when the palantiri_create command is called.


The main event configuration file is separated in several sub sections.

Clusterparameter
^^^^^^^^^^^^^^^^
The first subsection, [clusterparameter],  defines the parameters describing the cluster algorithm::

  [clusterparameter]

  The following parameter described the maximum number of desired array clusters that should be optimized for.
  Note that the actual number will vary from this.
  maxCluster = 100

  The number of the minimum number of stations in initial clusters (should be higher than the desired minimum number of stations per array)
  minStationAroundInitialCluster = 10

  The distance which is used initially to gather stations in an array. The final minimum array aperture will not depend on this. It should be  approximately 2x times the desired minimum array aperture.
  initialstationdistance = 7

  Cutoff is a hard break option in number of iterations for the clustering algorithm.
  cutoff = 100
  runs defines the number of times the clustering algorithm is applied.
  The final cluster setup to be used will be the one that performs best in terms of maximized number of stations and clusters used from all runs.
  runs = 1

  This defines the minimum distance of initial seach centroids [in degree]. Larger numbers means that less search centroids are spawned. This results in more distance between arrays.
  centroidminDistance = 2

  The maximum distance from station to cluster center (maximum array aperture)
  stationdistance = 5
  The minimum number of stations per cluster
  minClusterStation = 5


traveltime calculation options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The next chapter in the configuration file is for traveltime definitions and traveltime grid calculations::

  [traveltime calculation options]

  traveltime_model = ak135-f-continental.m.nd
  #number of cores for traveltime calculation
  ncore = 2

  #depths= Repeat semblance calculation from,to,steps relative to depth in origin config
  depths=0,0,5
  # run each depth step with a number of filter(s), used for high vs. low freq.
  filters=2

  # dimx of grid
  dimx = 20
  # dimy of grid
  dimy = 20

  True or false. The traveltimes and the semblance will be calculated in a cube. The parameter for the minimum and maximum depth and step size are given in the depths parameter.
  dimz= 0

  # gridspacing in degree
  gridspacing = 0.035
  # Phase to consider [right now only P possible!]
  ttphases=P

  futterman_attenuation = 0

Traveltime corrections
^^^^^^^^^^^^^^^^^^^^^^

traveltime correction options::

  [Traveltime corrections]

  Whether to apply shifts or not. Can be independent from the run to calculate the shifts.
  correct_shifts = 0

  Set true if you wish to calculate the empirical shifts or if you want to load in already calculated shifts.
  correct_shifts_empirical_run = 0

  Use empirical shift from fore or aftershock or calculate shifts based on maximizing the correlation [default].
  correct_shifts_empirical = 0

  Load in a manual file
  correct_shifts_empirical_manual = 0

  Calculate the empirical shifts per array (faster, default) or for each station individually.
  correct_shifts_empirical_manual_station_wise = 0

  Use synthetic or real event [default]:
  correct_shifts_empirical_synthetic=0

  # dimx of grid for empirical correction (should be the same as for the main process in most cases)
  dimx_emp = 1
  # dimy of grid (should be the same as for the main process in most cases)
  dimy_emp = 1
  # step length in s.
  step_emp = 16
  # window length in s.
  winlen_emp = 16
  the duration of the waveform in [s] after the theoretical onset of the specified phase
  duration_emp = 8
  the duration of the waveform in [s] before the theoretical onset of the specified phase
  forerun_emp = 8


Beamforming method
^^^^^^^^^^^^^^^^^^^^^^

Several methods for the backprojection are possible: 1. traditional (capon) beamforming
2. coherence based stacking (bp_coh), 3. phase weighted stacking (pws), 4. stacking in the frequency domain (bp_freq),
5. music (bp_music).
Compressed sensing is only partly supported right now.

Select the method for Backprojection in the following section::

  [beamforming method]

  #delaysum
  #capon
  beam = delaysum

  bp_freq= 0
  bp_coh= 0
  bp_music = 0

  # use a phase weighted stacking
  shift_by_phase_pws = 0

  # create output of compressed sensing as grid [warning: experimental]
  cs = 0

Algorithm settings
^^^^^^^^^^^^^^^^^^^^^
Define general settings::

  [algorithm settings]

  weight_by_noise = 0
  weight_by_azimuth = 0
  # bootstrap the arrays to estimate the uncertainty:
  bootstrap_array_weights = 0
  # number of bootstraps to carry out:
  n_bootstrap = 0

  Use arrays separately [default] or all in one big step:
  combine_all = 0
  Normalize the semblance of each array to 1.
  norm_all=1

General settings
^^^^^^^^^^^^^^^^^^^^^
Define general settings::

  [general parameter]

  # min distance to origin of stations
  minDist = 2
  # max distance to origin of stations
  maxDist = 93
  # if download of was done with palantiri_down command, set to 1
  pyrocko_download = 1
  # if download with pyrocko was done you can choose between velocity and displacement
  quantity = velocity
  Calculate and plot the array response of each array.
  array_response = 0

  Visualize the semblance after the calculation for direct inspection
  inspect_semb = 0


  [Synthetic Test]

  # do a synthetic test with a real station distribution, specify the
  # parameters in eventfolder with event.syn
  synthetic_test = 1
  # add noise to the synthetic, based on the variance of the real station
  # covariance of noise not enabled right now
  synthetic_test_add_noise = 0
  synthetic_test_pertub_arrivals = 0
  shift_max = 10
  # if colosseo synthetics should be used, set to 1
  colesseo_input = 0
  # give the colosseo scenario.yml file
  colosseo_scenario_yml = /media/asteinbe/data/asteinbe/mydl/scenario.yml


  [Filter and time window settings]

  # step length in s.
  step = 6
  # window length in s.
  winlen = 12
  # step length in s.
  step_f2 = 10
  # window length in s.
  winlen_f2 = 20
  # length of data before phase onset in s.
  forerun = 10
  # length of data after phase onset in s.
  duration = 40
  # resampling data to frequency in Hz or s, should match your gf store
  new_frequence = 0.5


  [Manual shifting]

  # shift the traces to theoretical onset
  shift_by_phase_onset = 0
  # shift by crosscorrelation
  shift_by_phase_cc = 0


  [Optimization]


  # Optimize for pyrocko sources with array responses as input with the semblance(all)
  optimize = 0
  optimize_all = 0

  [focal mechanism solution values from event file]
  #only = 1 possible
  fm = 1

  [xcorrskript parameter]

  xcorr=0
  # for manual qc set autoxcorrcorrectur to 1:
  autoxcorrcorrectur = 0
  # crosscorrelation threshold for excluding stations
  xcorrtreshold = 0.6

  #filter for referencestation for automatic picker
  #should match your filter
  refstationfreqmin=0.03
  refstationfreqmax=0.08
  refstationcorners=2
  refstationzph=false

  #STA/LTA parameter
  refsta=0.5
  reflta=4


  [filterparameter]

  filterswitch=1
  Calculate the filters based on the estimated corner frequency. The source parameters need to be defined in the event config file.
  dynamic_filter = 0

  ###############################################################
  #Parameter for first filter
  #bp butterworth

  # low cut corner frequency
  flo = 0.1

  # high cut corner frequency
  fhi = 0.24

  # number of filter sections
  ns = 4

  # TRUE -> zero phase filter
  zph = false


  ###############################################################
  #Example Parameter for a second filter
  #bp butterworth

  # low cut corner frequency
  flo2 = 0.03

  # high cut corner frequency
  fhi2 = 0.24

  # number of filter sections
  ns2 = 4

  # TRUE -> zero phase filter
  zph2 = false


  [array parameter]
  Here follow all the arrays. This will be updated by the cluster algorithm.

  networks=r1,r2
  r1=XK.B03SL..BHZ|XK.B04KH..BHZ|XK.B05MO..BHZ|XK.B06OR..BHZ|XK.B07DX..BHZ|XK.B08TS..BHZ|XK.B09NK..BHZ|XK.B10PP..BHZ|XK.B11ET..BHZ|XK.B12SS..BHZ|XK.B13NX..BHZ|XK.B14MH..BHZ|XK.B15MW..BHZ|XK.B17CI..BHZ
  r1refstation=
  r1phase=P
  r2=NM.OLIL..BHZ|PN.PPCWF..BHZ|TA.O44A..BHZ|TA.O45A..BHZ|TA.P43A..BHZ|TA.P44A..BHZ|TA.P45A..BHZ|TA.P46A..BHZ|TA.Q43A..BHZ|TA.Q44A..BHZ|TA.Q45A..BHZ|TA.Q46A..BHZ|TA.R45A..BHZ|XO.LA19..BHZ|XO.LA21..BHZ|XO.LB20..BHZ|XO.LB22..BHZ
  r2refstation=
  r2phase=P
