import logging
import fnmatch
import os
import shutil
import glob
import numpy as np
import sys
from pyrocko import model
from pyrocko.guts import Object, Float, Int, String, Bool, List, Tuple

if sys.version_info.major >= 3:
    from configparser import SafeConfigParser
else:
    from ConfigParser import SafeConfigParser

logger = logging.getLogger('ARRAY-MP')


class PalantiriDataConfig(Object):
    '''Configuration of data IO and data preprocessing'''
    quantity = String.T(default="displacement",
                            help='quantity')

    tttopt = String.T(default="-ph P",
                            help='quantity')

    cam = Int.T(default=2, help='Length in seconds. Not needed \
        when using TFRecordData')

    pyrocko_download = Bool.T(default=True, optional=True,
     help='Length in seconds. Not needed \
        when using TFRecordData')

    export_unfiltered = Bool.T(default=False, optional=True,
     help='')

    export_filtered = Bool.T(default=False, optional=True,
     help='')

    export_resampled = Bool.T(default=False, optional=True,
     help='')

    colesseo_input = Bool.T(default=False, optional=True,
     help='')

    colosseo_scenario_yml = String.T(optional=True, help='')

    load_wdf = Bool.T(default=False, optional=True,
     help='')


class PalantiriSyntheticConfig(Object):
    '''Configuration of data IO and data preprocessing'''

    synthetic_test = Bool.T(default=False, optional=True,
     help='')

    synthetic_test_add_noise = Bool.T(default=False, optional=True,
     help='')

    synthetic_test_pertub_arrivals = Bool.T(default=False, optional=True,
     help='')


class PalantiriWeightConfig(Object):
    '''Configuration of data IO and data preprocessing'''

    shift_max = Float.T(default=4., optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    weight_by_azimuth = Bool.T(default=True, optional=True,
     help='')

    bootstrap_array_weights = Bool.T(default=True, optional=True,
     help='')

    n_bootstrap = Int.T(default=1, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    correct_shifts_empirical_run = Bool.T(default=False, optional=True,
     help='')

    correct_shifts = Bool.T(default=False, optional=True,
     help='')

    correct_shifts_empirical = Bool.T(default=False, optional=True,
     help='')

    correct_shifts_empirical_manual = Bool.T(default=False, optional=True,
     help='')

    correct_shifts_empirical_manual_station_wise = Bool.T(default=False, optional=True,
     help='')

    combine_all = Bool.T(default=False, optional=True,
     help='')

    norm_all = Bool.T(default=True, optional=True,
     help='')

    weight_by_noise = Bool.T(default=False, optional=True,
     help='')

    shift_by_phase_onset = Bool.T(default=False, optional=True,
     help='')

    shift_by_phase_pws = Bool.T(default=False, optional=True,
     help='')

    shift_by_phase_cc = Bool.T(default=False, optional=True,
     help='')


class PalantiriFilterConfig(Object):
    '''Configuration of data IO and data preprocessing'''

    filters = Int.T(default=1, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    dynamic_filter = Bool.T(default=False, optional=True,
     help='')

    filterswitch = Int.T(default=1, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    flo = List.T(
        Float.T(), help='List blacklist patterns (may contain wild cards')

    fhi = List.T(
        Float.T(), help='List blacklist patterns (may contain wild cards')

    name = List.T(
        String.T(), help='List blacklist patterns (may contain wild cards')

    #zph = List.T(
    #    String.T(), help='List blacklist patterns (may contain wild cards')

    newFrequency = Float.T(default=0.5, help='Length in seconds. Not needed \
        when using TFRecordData')

    ns = List.T(
        Int.T(), help='List blacklist patterns (may contain wild cards')

    step_emp = Float.T(default=4., optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    winlen_emp = Float.T(default=4., optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    step = Float.T(default=2., optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    winlen = Float.T(default=8., optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    forerun = Float.T(default=10., optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    duration = Float.T(default=20., optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')


class PalantiriGeometryConfig(Object):
    '''Configuration of data IO and data preprocessing'''

    dimy = Int.T(default=50, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    dimx = Int.T(default=50, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    dimz = Int.T(default=0, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    gridspacing = Float.T(default=0.025, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    dimx_emp = Int.T(default=50, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    dimy_emp = Int.T(default=50, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    depths = List.T(
        Float.T(), help='List blacklist patterns (may contain wild cards')


class PalantiriXcorrConfig(Object):
    '''Configuration of data IO and data preprocessing'''

    xcorr = Bool.T(default=True, optional=True,
                   help='')

    xcorrtreshold = Float.T(default=0.6, optional=True, help='')

    autoxcorrcorrectur = Bool.T(default=False, optional=True,
                                help='')

    refstationfreqmin = Float.T(default=0.03, optional=True,
                                help='Length in seconds. Not needed \
                                when using TFRecordData')

    refstationfreqmax = Float.T(default=1, optional=True,
                                help='Length in seconds. Not needed \
                                when using TFRecordData')

    refstationcorners = Int.T(default=2, optional=True,
                              help='Length in seconds. Not needed \
                              when using TFRecordData')

    refstationzph = Bool.T(default=False, optional=True, help='')

    refsta = Float.T(default=0.5, optional=True,
                     help='Length in seconds. Not needed \
                     when using TFRecordData')

    reflta = Float.T(default=2, optional=True,
                     help='Length in seconds. Not needed \
                     when using TFRecordData')


class PalantiriConfig(Object):
    '''Configuration of data IO and data preprocessing'''

    blacklist = List.T(
        String.T(), help='List blacklist patterns (may contain wild cards')

    inspect_semb = Bool.T(default=False, optional=True,
     help='')

    traveltime_model = String.T(optional=True, default="ak135-f-continental.m.nd", help='')

    futterman_attenuation = Bool.T(default=False, optional=True,
     help='')

    bp_freq = Bool.T(default=False, optional=True,
     help='')

    bp_coh = Bool.T(default=False, optional=True,
     help='')

    bp_music = Bool.T(default=False, optional=True,
     help='')

    optimize = Bool.T(default=False, optional=True,
     help='')

    optimize_all = Bool.T(default=False, optional=True,
     help='')

    cs = Bool.T(default=False, optional=True,
     help='')

    array_response = Bool.T(default=False, optional=True,
     help='')

    ncores = Int.T(default=2, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    beam = String.T(default="delaysum", optional=True,
                            help='quantity')

    ttphases = List.T(
        String.T(), help='List blacklist patterns (may contain wild cards')


class ClusterConfig(Object):
    '''Configuration of data IO and data preprocessing'''

    cluster = List.T(
        String.T(), help='List blacklist patterns (may contain wild cards')

    maxcluster = Int.T(default=100, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    minStationAroundInitialcluster = Int.T(default=4, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    cutoff = Int.T(default=10, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    runs = Int.T(default=1, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    centroidminDistance = Float.T(default=2, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    comparedelta = Float.T(default=2, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    stationdistance = Float.T(default=5, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    minclusterStation = Float.T(default=4, optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    minDist = Float.T(default=23., optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    initialstationdistance = Float.T(default=5., optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')

    maxDist = Float.T(default=93., optional=True, help='Length in seconds. Not needed \
        when using TFRecordData')


class Config(Object):
    '''Defines your setup.'''

    name = String.T(default='unnamed',
        help='Used to identify the model and runs in summy and checkpoint \
            directory')

    config = PalantiriConfig.T(help='')
    config_data = PalantiriDataConfig.T(help='')
    config_syn = PalantiriSyntheticConfig.T(help='')
    config_weight = PalantiriWeightConfig.T(help='')
    config_filter = PalantiriFilterConfig.T(help='')
    config_geometry = PalantiriGeometryConfig.T(help='')
    config_xcorr = PalantiriXcorrConfig.T(help='')
    config_cluster = ClusterConfig.T(help='')


class Station(object):
    '''
    class to store station object from metadatafile
    '''

    def __init__(self, net, sta, loc, comp, lat=0, lon=0, ele=0, dip=0, azi=0,
                 gain=0, takeoff=0, backazi=0, sazi=0):
        self.net = net
        self.sta = sta
        self.loc = loc
        self.comp = comp
        self.lat = lat
        self.lon = lon
        self.ele = ele
        self.dip = dip
        self.azi = azi
        self.gain = gain
        self.takeoff = takeoff
        self.backazi = backazi
        self.sazi = sazi

    def getName(self):
        return (self.net+'.'+self.sta+'.'+self.loc+'.'+self.comp)

    def pyr_name(self):
        return self.net+'.'+self.sta+'.'+self.loc+'.'

    def getcmpName(self):
        if self.loc == '--':
            self.loc = ''
        return self.net+'.'+self.sta+'.'+self.loc+'.'+self.comp

    def __str__(self):
        return('%s.%s.%s.%s') % (self.net, self.sta, self.loc, self.comp)

    def __eq__(self, other):
        return self.getName() == other.getName()


class Event(object):
    '''
    class to store event object form origin file
    '''
    def __init__(self, lat, lon, depth, time='', region='', strike=0,
                 dip=0, rake=0):
        self.lat = lat
        self.lon = lon
        self.depth = depth
        self.region = region
        self.time = time
        self.strike = strike
        self.dip = dip
        self.rake = rake


class Trigger(object):
    '''
    class to store triggertimes for station used in xcorrelations
    '''

    def __init__(self, stationname, triggertime, arrayname, tdiff=0):
        self.sname = stationname
        self.ttime = triggertime
        self.aname = arrayname
        self.tdiff = tdiff


class Config(object):
    '''
    class to parse origin, config and metadata file for arraytool
    initialize with eventpath
    '''

    def __init__(self, eventpath, eventpath_emp=None):
        self.eventpath = eventpath
        self.eventpath_emp = eventpath_emp

    def parseConfig(self, suffix):
        '''
        method to parse config files(origin,config)
        return Configdictionary
        '''
        cDict = {}
        if suffix == "origin_emp":
            try:
                files = glob.glob(os.path.join(self.eventpath_emp, '*.'+'origin'))
                parser = SafeConfigParser()
                parser.read(files[0])
            except:
                files = glob.glob(os.path.join(self.eventpath_emp, '*.'+'origin_emp'))
                parser = SafeConfigParser()
                parser.read(files[0])
            for section_name in parser.sections():
                for name, value in parser.items(section_name):
                    cDict[name] = value
        elif suffix == "yaml" or suffix == "yml":
            try:
                cDict = glob.glob(os.path.join(self.eventpath, '*.'+'yaml'))
            except:
                cDict = glob.glob(os.path.join(self.eventpath, '*.'+'yml'))
        else:
            files = glob.glob(os.path.join(self.eventpath, '*.'+suffix))
            parser = SafeConfigParser()
            parser.read(files[0])

            for section_name in parser.sections():
                for name, value in parser.items(section_name):
                    cDict[name] = value

        return cDict

    def readMetaInfoFile(self):
        '''
        method to parse metadata file
        return List of Station Objects
        '''
        MetaL = []
        logger.info('\033[31m Parsing MetaInfoFile \033[0m \n')
        try:
            for i in os.listdir(self.eventpath):
                if fnmatch.fnmatch(i, '*.meta'):
                    evfile = os.path.join(self.eventpath, i)
                    fobj = open(evfile, 'r')
                    for i in fobj:
                        line = i.split()
                        net = line[0]
                        sta = line[1]
                        loc = line[2]
                        comp = line[3]
                        lat = line[4]
                        lon = line[5]
                        ele = line[6]
                        dip = line[7]
                        azi = line[8]
                        gain = line[9]
                        if fnmatch.fnmatch(comp, '*HZ'):
                            MetaL.append(Station(net, sta, loc, comp, lat, lon,
                                                 ele, dip, azi, gain))

            logger.info('\033[31m %d ENTRIES IN METAFILE FOUND \033[0m \n' % (len(MetaL)))
        except:
            logger.info('\033[31m METAFILE NOT READABLE \033[0m \n')

        FML = self.checkMetaInfoFile(MetaL)

        return FML

    def readpyrockostations(self, phases, test=False):
        try:
            stations = model.load_stations(self.eventpath+'/data/stations_cluster.txt')
        except Exception:
            stations = model.load_stations(self.eventpath+'/data/stations_disp.txt')
        MetaL = []
        for phase in phases:
            if phase is 'P':
                desired = 'Z'
            if phase is 'S':
                desired = 'T'
            for sl in stations:
                count_channel = 0
                for channel in sl.channels:
                    if channel.name[-1] == desired and channel.name is not "HHZ":
                        MetaL.append(Station(str(sl.network), str(sl.station),
                                             str(sl.location),
                                             str(channel)[:3], str(sl.lat),
                                             str(sl.lon),
                                             str(sl.elevation),
                                             str(channel.dip),
                                             str(channel.azimuth),
                                             str(channel.gain)))
                    count_channel = count_channel+1

        FML = self.checkMetaInfoFile(MetaL)

        return FML

    def readcolosseostations(self, scenario_path):
        stations = model.load_stations(scenario_path+'/meta/stations.txt')

        MetaL = []
        for sl in stations:
                channel = sl.channels[2]
                MetaL.append(Station(str(sl.network), str(sl.station),
                             str(sl.location), str(sl.channels[2])[:3],
                             str(sl.lat), str(sl.lon), str(sl.elevation),
                             str(channel.dip), str(channel.azimuth),
                             str(channel.gain)))
        FML = self.checkMetaInfoFile(MetaL)

        return FML

    def checkMetaInfoFile(self, MetaList):

        ML = []
        DL = []
        LL = []

        for i in MetaList:
            try:
                if float(i.gain) == 0:
                    print('GAIN IS ZERO ', i)
                    search = ('%s.%s.%s') % (i.net, i.sta, i.loc)
                    DL.append(search)
                    LL.append(i)
            except:
                i.gain = (np.float(i.gain[:-3]))  # careful, there is something off with some japanese/chinese stats.
                search = ('%s.%s.%s') % (i.net, i.sta, i.loc)
                DL.append(search)
                LL.append(i)

        if len(DL) > 0:
            for i in DL:
                for j in MetaList:
                    metaname = ('%s.%s.%s') % (j.net, j.sta, j.loc)
                    if i != metaname:
                        ML.append(j)
        else:
            ML = MetaList
        return ML

    def createFolder(self):
        '''
        method to create work folder in event directory
        returns Folderdictionary
        '''
        Folder = {}
        logger.info('\033[31m Create working environment \033[0m \n')
        if os.access(os.getcwd(), os.W_OK):

            basedir = os.path.join(self.eventpath, 'work')
            sembdir = os.path.join(basedir, 'semblance')
            ascdir = os.path.join(basedir, 'asc')
            mseeddir = os.path.join(basedir, 'mseed')

            Folder['base'] = basedir
            Folder['semb'] = sembdir
            Folder['asc'] = ascdir
            Folder['mseed'] = mseeddir

            for key in Folder:
                if os.access(Folder[key], os.F_OK) is False:
                    os.makedirs(Folder[key])

        else:
            print("no write permissions")

        Folder['config'] = os.path.join('..', 'skeleton')
        return Folder

    def writeConfig(self, Config, Origin, Folder):
        '''
        method to write recently used config to event folder
        '''
        dst_dir = os.path.join(Folder['semb'], 'config_file.cfg')

        files = glob.glob(os.path.join(self.eventpath, '*.'+'config'))
        parser = SafeConfigParser()
        src_file = parser.read(files[0])

        shutil.copy(src_file[0], dst_dir)

    def writeStationFile(self, StationMetaData, Folder, flag):
        '''
        method to write recently used stations for processing to event folder
        '''
        name = 'stations_'+str(flag)+'.dat'
        fobj = open(os.path.join(Folder['semb'], name), 'w')
        for i in StationMetaData:
            fobj.write('%s %s %s %s\n' % (flag, i.getName(), i.lat, i.lon))
        fobj.close()
