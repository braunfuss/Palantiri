#
#    Common data classes
#
# -------------------------------------------------------------------------------------------------

class Station(object):

    def __init__(self, net, sta, loc, comp,lat=0,lon=0,ele=0,dip=0,azi=0,gain=0,inst=None):

        self.net  = net
        self.sta  = sta
        self.loc  = loc
        self.comp = comp
        self.site = '???'          # fuer getStationList.py  ???
        self.lat  = lat
        self.lon  = lon
        self.ele  = ele
        self.dip  = dip
        self.azi  = azi
        self.gain = gain
        self.inst = inst
        self.provider = None

    def fullName(self) :    return self.net + '.' + self.sta + '.' + self.loc + '.' + self.comp  
    def stationName(self) : return self.net + '.' + self.sta
    def location(self) :    return Location(self.lat, self.lon)

    def print1(self) :

        print 'name = ', self.fullName()
        print 'lat  = ', self.lat
        print 'lon  = ', self.lon
        print 'ele  = ', self.ele
        print 'dip  = ', self.dip
        print 'azi  = ', self.azi
        print 'gain = ', self.gain
        print 'inst = ', self.inst

        if self.provider != None : print 'prov = ' + self.provider

        print ' '
#endClass

# -------------------------------------------------------------------------------------------------

class Location(object) :

    def __init__(self, lat, lon):
 
        self.lat  = float(lat)
        self.lon  = float(lon)

    def __str__(self):        return('(%f,%f)') %(self.lat,self.lon)
    def __eq__(self, other) : return(self.lat == other.lat and self.lon == other.lon)

    def set(d) :

        self.lat = d.lat
        self.lon = d.lon

#endClass

def dictToLocation(d) :

    lat = float(d['lat'])
    lon = float(d['lon'])
    return Location(lat, lon)

# -------------------------------------------------------------------------------------------------

def _getDelim(name) :

    if   '.' in name : return '.'
    elif '_' in name : return '_'

    assert False


def toStationNames(strings) :

    names = []

    if len(strings) == 0 :  return names
    
    delim = _getDelim(strings[0])

    for i in range(len(strings)) :
       s = strings[i].split(delim)
       names.append(s[0] + delim + s[1])

    return names
# -------------------------------------------------------------------------------------------------

def toNetwork(name) :
    s = name.split(_getDelim(name))
    return s[0]

def toStation(name) :
    s = name.split(_getDelim(name))
    return s[1]

def toNetAndStation(name) :
    return toNetwork(name) + _getDelim(name) + toStation(name)

def isSameNetwork(s1,s2) :
    return(toNetwork(s1) == toNetwork(s2))

def toNetworkNames(strings) :

    names = []

    for i in range(len(strings)) :
       names.append(toNetwork(strings[i]))

    return names

def selectNetwork(stationList, network) :

    result = []

    for s in stationList :
        if s.startswith(network) : result.append(s)

    return result

# -------------------------------------------------------------------------------------------------
