from configparser import SafeConfigParser
from urllib.parse import urlencode
from urllib.request import urlopen


class Event(object):
    def __init__(self, id, region, lat, lon, otime, mag):
        self.id = id
        self.region = region
        self.lat = lat
        self.lon = lon
        self.otime = otime
        self.mag = mag

    def __str__(self):
        return 'EventID: {0:10} ---> {1:35} {2:30} {3:10}{4:12}{5:12}'.format(
            self.id, self.region, self.otime, self.mag, self.lat, self.lon)


def init():
    cDict = {}
    parser = SafeConfigParser()
    parser.read('global.conf')
    for section_name in parser.sections():
        for name, value in parser.items(section_name):
            cDict[name] = value
    return cDict


def parseIrisEventWebservice(searchparameter):

    if not searchparameter['resultlimit']:
        searchparameter['resultlimit'] = 10

    url = 'http://service.iris.edu/fdsnws/event/1/query?'
    catalog = searchparameter['catalog']
    if len(catalog) > 1:
            parameter = urlencode({
                         'catalog': searchparameter['catalog'],
                         'minmag': searchparameter['magmin'],
                         'starttime': searchparameter['date_min'],
                         'endtime': searchparameter['date_max'],
                         'format': 'text',
                         'orderby': 'time',
                         'limit': searchparameter['resultlimit'],
                    })
    else:
        parameter = urlencode({
                     'minmag': searchparameter['magmin'],
                     'starttime': searchparameter['date_min'],
                     'endtime': searchparameter['date_max'],
                     'format': 'text',
                     'orderby': 'time',
                     'limit': searchparameter['resultlimit'],
                })
    u = ('%s%s') % (url, parameter)

    data = urlopen(u).read()
    data = data.decode('utf_8')
    data = data.split('\n')
    dl = data[1:]

    EL = []
    for i in dl:
        if len(i) != 0:
            i = i.split('|')
            EL.append(Event(i[0], i[12], i[2], i[3], i[1], i[10]))

    print('\n\n         # to get data for event use eventID #\n\n')
    if len(EL) is not 0:
        for event in EL:
            print(event)
    else:
        print('\033[31m No event entry found \033[0m\n')


def searchEvent(searchparameter):

    parseIrisEventWebservice(searchparameter)


def main():
    options = init()
    searchEvent(options)
