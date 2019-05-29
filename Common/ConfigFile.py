import os
import sys
import Globals
import Basic
import Logfile
if sys.version_info.major >= 3:
    from configparser import SafeConfigParser
else:
    from ConfigParser import SafeConfigParser


class ConfigObj(object):

    def __init__(self, fileName=None, dict=None):
        self._fileName = fileName
        self._dict = dict
        self._debug = False

    def __getitem__(self, key):
        return self._dict.__getitem__(key)

    def getDict(self): return self._dict

    def setDict(self, d): self._dict = d

    def getFileName(self): return self._fileName

    def Str(self, key, default=None):
        return self.String(key, default)

    def String(self, key, default=None):

        if key in self._dict:
            s = self._dict[key]
        elif default is not None:
            s = default
        else:
            self._keyMissing(key)

        if self._debug:
            print('cfg [' + key + '] = ' + s)

        return s

    def Int(self, key, default=None):

        if default is not None:
            val = self._checkIntNumber(key, default)
        else:
            val = None

        s = self.String(key, val)
        return self._checkIntNumber(key, s)

    def UInt(self, key, maxVal=None, default=None):

        assert default is None or default >= 0
        val = self.Int(key, default)
        return val

    def Float(self, key, default=None):

        if default is not None:
            val = self._checkFloatNumber(key, default)
        else:
            val = None

        s = self.String(key, val)
        return self._checkFloatNumber(key, s)


    def UFloat(self, key, maxVal=None, default=None):

        assert default is None or default >= 0.0
        val = self.Float(key, default)

        if val < 0.0:
            self._error(key, str(val) + ' < 0.0')
        if maxVal is not None and val > maxVal:
            self._error(key, str(val) + ' > ' + str(maxVal))

        return val

    def Bool(self, key, default=None):

        val = self.Int(key, default)

        if val == 0:
            return False
        if val == 1:
            return True

        self.rangeError(key, '0', '1')

    def IntRange(self, key1, key2):

        a = self.Int(key1)
        b = self.Int(key2)

        if a > b:
            self.rangeError2(key1, key2)

        return a, b

    def FloatRange(self, key1, key2):

        a = self.Float(key1)
        b = self.Float(key2)

        if a > b:
            self.rangeError2(key1, key2)

        return a, b


    def Distance(self, key, maxVal=40 * 1000):
        return self.UFloat(key, maxVal)

    def Freq(self, key, maxVal=1000):
        return self.UFloat(key, maxVal)

    def Duration(self, key='duration'):
        return self.UInt(key)

    def Time(self, key='time'):
        return self.String(key)

    def traveltime_model(self):
        return self.Str('traveltime_model')

    def lat(self):
        return self.Float('lat')

    def lon(self):
        return self.Float('lon')

    def depth(self):
        return self.Float('depth')

    def duration_emp(self):
        return self.Float('duration_emp')

    def forerun_emp(self):
        return self.Float('forerun_emp')

    def inspect_semb(self):
        return self.Bool('inspect_semb')

    def array_response(self):
        return self.Bool('array_response')

    def dimX(self):
        return self.UInt('dimx')

    def dimY(self):
        return self.UInt('dimy')

    def dimX_emp(self):
        return self.UInt('dimx_emp')

    def dimY_emp(self):
        return self.UInt('dimy_emp')

    def winlen_emp(self):
        return self.Float('winlen_emp')

    def step_emp(self):
        return self.Float('step_emp')

    def winlen(self):
        return self.Float('winlen')

    def step(self):
        return self.Float('step')

    def winlen_f2(self):
        return self.Float('winlen_f2')

    def step_f2(self):
        return self.Float('step_f2')

    def cs(self):
        return self.Bool('compressed_sensing')

    def synthetic_test(self):
        return self.Bool('synthetic_test')

    def shift_by_phase_onset(self):
        return self.Bool('shift_by_phase_onset')

    def shift_by_phase_pws(self):
        return self.Bool('shift_by_phase_pws')

    def shift_by_phase_cc(self):
        return self.Bool('shift_by_phase_cc')

    def newFrequency(self):
        return self.Float('new_frequence')

    def pyrocko_download(self):
        return self.Bool('pyrocko_download')

    def synthetic_test_add_noise(self):
        return self.Bool('synthetic_test_add_noise')

    def synthetic_test_pertub_arrivals(self):
        return self.Bool('synthetic_test_pertub_arrivals')

    def shift_max(self):
        return self.Float('shift_max')

    def quantity(self):
        return self.Str('quantity')

    def ttphases(self):
        return self.Str('ttphases')

    def colesseo_input(self):
        return self.Bool('colesseo_input')

    def optimize(self):
        return self.Bool('optimize')

    def optimize_all(self):
        return self.Bool('optimize_all')

    def combine_all(self):
        return self.Bool('combine_all')

    def colosseo_scenario_yml(self):
        return self.Str('colosseo_scenario_yml')

    def weight_by_noise(self):
        return self.Bool('weight_by_noise')

    def norm_all(self):
        return self.Bool('norm_all')

    def futterman_attenuation(self):
        return self.Bool('futterman_attenuation')

    def dynamic_filter(self):
        return self.Bool('dynamic_filter')

    def bootstrap_array_weights(self):
        return self.Bool('bootstrap_array_weights')

    def n_bootstrap(self):
        return self.UInt('n_bootstrap')

    def weight_by_azimuth(self):
        return self.Bool('weight_by_azimuth')

    def correct_shifts(self):
        return self.Bool('correct_shifts')

    def correct_shifts_empirical(self):
        return self.Bool('correct_shifts_empirical')

    def correct_shifts_empirical_run(self):
        return self.Bool('correct_shifts_empirical_run')


    def _error0(self, msg):

        Logfile.error(msg)
        Logfile.abort()

    def _error(self, key, msg):
        self._error0(self.getFileName() + ',' + key + ': ' + msg)

    def _checkFloatNumber(self, key, s):

        if not Basic.isNumber(s):
            self._error(key, 'Key is not a number')
        return float(s)

    def _checkIntNumber(self, key, s):

        if not Basic.isNumber(s):
            self._error(key, 'Key is not a number')
        if not Basic.isInt(s):
            self._error(key, 'Key is not integer number')

        return int(s)

    def _keyMissing(self, key):
        self._error(key, 'Key missing')

    def rangeError(self, key, a, b):
        self._error(key, 'Value outside [' + a + ',' + b + ']')

    def rangeError2(self, key1, key2):
        self._error0('Range error: ' + key1 + ' > ' + key2)


class FilterCfg(ConfigObj):

    def __init__(self, dict):   ConfigObj.__init__(self, None, dict)

    def newFrequency(self):
        return self.Float('new_frequence')

    def filterswitch(self):
        return self.UInt('filterswitch', 3)

    def flo(self):
        return self.Float('flo')

    def fhi(self):
        return self.Float('fhi')

    def ns(self):
        return self.Float('ns')

    def flo2(self):
        return self.Float('flo2')

    def fhi2(self):
        return self.Float('fhi2')

    def ns2(self):
        return self.Float('ns2')

    def l_fc(self):
        return self.Freq('l_fc')

    def l_ns(self):
        return self.Freq('l_ns')

    def h_fc(self):
        return self.Freq('h_fc')

    def h_ns(self):
        return self.Freq('h_ns')

    def filterName(self):

        i = self.filterswitch()

        if i == 0:
            return None
        if i == 1:
            strings = [str(self.flo()),  str(self.fhi()),
                       str(self.ns()), self.Str('zph')]
        elif i == 2:
            strings = [str(self.l_fc()), str(self.l_ns()), self.Str('l_zph')]
        elif i == 3:
            strings = [str(self.h_fc()), str(self.h_ns()), self.Str('h_zph')]

        name = '_'.join(strings)
        return name


def filterName(dict):

    cfg = FilterCfg(dict)
    return cfg.filterName()


class OriginCfg(ConfigObj):

    def __init__(self, dict):
        ConfigObj.__init__(self, None, dict)

    def strike(self, def1):
        return self.Float('strike', def1)

    def dip(self, def1):
        return self.Float('dip',    def1)

    def rake(self, def1):
        return self.Float('rake',   def1)

    def time(self):
        return self.String('time')


class SynthCfg(ConfigObj):

    def __init__(self, dict):   ConfigObj.__init__(self, None, dict)

    def lat_0(self): return self.Float('lat_0')
    def lon_0(self): return self.Float('lon_0')
    def north_shift_0(self): return self.Float('north_shift_0')
    def east_shift_0(self): return self.Float('east_shift_0')
    def strike_0(self): return self.Float('strike_0')
    def dip_0(self): return self.Float('dip_0')
    def rake_0(self): return self.Float('rake_0')
    def length_0(self): return self.Float('length_0')
    def width_0(self): return self.Float('width_0')
    def depth_syn_0(self): return self.Float('depth_0')
    def nucleation_x_0(self): return self.Float('nucleation_x_0')
    def nucleation_y_0(self): return self.Float('nucleation_y_0')
    def slip_0(self): return self.Float('slip_0')
    def rmnn(self): return self.Float('rmnn')
    def rmee(self): return self.Float('rmee')
    def rmdd(self): return self.Float('rmdd')
    def rmne(self): return self.Float('rmne')
    def rmnd(self): return self.Float('rmnd')
    def rmed(self): return self.Float('rmed')
    def duration(self): return self.Float('duration')
    def magnitude_0(self): return self.Float('magnitude_0')
    def mag_0_low(self): return self.Float('mag_low_0')
    def mag_0_high(self): return self.Float('mag_high_0')
    def dip_0_low(self): return self.Float('dip_low_0')
    def dip_0_high(self): return self.Float('dip_high_0')
    def depth_0_low(self): return self.Float('depth_low_0')
    def depth_0_high(self): return self.Float('depth_high_0')
    def rake_0_low(self): return self.Float('rake_low_0')
    def rake_0_high(self): return self.Float('rake_high_0')
    def strike_0_low(self): return self.Float('strike_low_0')
    def strike_0_high(self): return self.Float('strike_high_0')
    def north_shift_0_low(self): return self.Float('north_shift_low_0')
    def north_shift_0_high(self): return self.Float('north_shift_high_0')
    def east_shift_0_low(self): return self.Float('east_shift_low_0')
    def east_shift_0_high(self): return self.Float('east_shift_high_0')
    def time_0_low(self): return self.Float('time_low_0')
    def time_0_high(self): return self.Float('time_high_0')

    def time_0(self): return self.String('time_0')
    def store_superdirs(self): return self.Str('store_superdirs')
    def store(self): return self.Str('store')
    def source(self): return self.Str('source')
    def stf(self): return self.Str('stf')
    def use_specific_stf(self): return self.Bool('use_specific_stf')
    def nsources(self): return self.Int('nsources')
    def stf(self): return self.Str('stf')

    def lat_1(self, i): return self.Float(('lat_%s' % i))
    def lon_1(self, i): return self.Float(('lon_%s' % i))
    def north_shift_1(self, i): return self.Float(('north_shift_%s' % i))
    def east_shift_1(self, i): return self.Float(('east_shift_%s' % i))
    def strike_1(self, i): return self.Float(('strike_%s' % i))
    def depth_syn_1(self, i): return self.Float(('depth_%s' % i))
    def dip_1(self, i): return self.Float(('dip_%s' % i))
    def rake_1(self, i): return self.Float(('rake_%s' % i))
    def length_1(self, i): return self.Float(('length_%s' % i))
    def width_1(self, i): return self.Float(('width_%s' % i))
    def nucleation_x_1(self, i): return self.Float(('nucleation_x_%s' % i))
    def nucleation_y_1(self, i): return self.Float(('nucleation_y_%s' % i))
    def slip_1(self, i): return self.Float(('slip_%s' % i))
    def magnitude_1(self, i): return self.Float(('magnitude_%s' % i))
    def time_1(self, i): return self.String(('time_%s' % i))


DEFAULT_CONFIG_FILE = 'global.conf'
blacklist = 'blacklist'
duration = 'duration'
keyfilefolder = 'keyfilefolder'

mail = 'mail'
mindist = 'mindist'
maxdist = 'maxdist'

metaCatalog = 'metacatalog'
pwd = 'pwd'


class GlobalConfigObj(ConfigObj):

    def __init__(self, fileName):

        if fileName is None:
            name = DEFAULT_CONFIG_FILE
        else:
            name = fileName

        ConfigObj.__init__(self, name, readConf(os.path.join('..', name)))


_globConfigObj = None

def GlobCfg():
    return _globConfigObj

def readGlobalConf(fileName):
    global _globConfigObj

    _globConfigObj = GlobalConfigObj(fileName)
    return _globConfigObj.getDict()

# -------------------------------------------------------------------------------------------------

def readConf(fileName):

    if not Basic.checkFileExists(fileName):
        return None

    cDict = {}
    parser = SafeConfigParser()
    parser.read(fileName)

    isClient = Globals.isClient

    if not isClient:
        Logfile.setVisible(False)
        Logfile.add(' ', fileName, ': ')

    Logfile.setErrorLog(True)

    for section_name in parser.sections():
        for name, value in parser.items(section_name):
            cDict[name] = value

            if not isClient:
                if name != 'mail' and name != 'pwd':
                    Logfile.add(name + ' = ' + value)

    if not isClient:
        Logfile.add(' ')
        Logfile.setVisible(True)

    Logfile.setErrorLog(False)
    return cDict


def checkKeys(conf, keyList, optional=False):

    if type(keyList) is str:
        list1 = list(keyList)
    else:
        list1 = keyList

    if not optional:
        Basic.checkExistsKeys(conf, list1, isAbort=True)

    eventDir = Globals.EventDir()
    isOk = True

    for key in list1:
        val = conf[key]
        msg = None

        if key == duration:
            msg = Basic.checkGreaterZero(val)
        elif key in [mindist, maxdist]:
            msg = Basic.checkNotNegative(val)
        elif key in [keyfilefolder, metaCatalog]:
            Basic.checkExistsDir(os.path.join(eventDir, val), isAbort=True)
        elif key in [blacklist, mail, pwd]:
            continue

        if msg is not None:
            isOk = Logfile.error('Key <' + key + '> in config file: ' + msg)

    if not isOk:
        Logfile.abort()
    return True
