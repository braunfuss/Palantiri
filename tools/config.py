import logging
import os.path as op
from optparse import OptionParser

from pyrocko import util, scenario, guts, gf
import fnmatch
import os
import shutil
import glob
import logging
import platform
import numpy as np
from ConfigParser import SafeConfigParser
from pyrocko import model
logger = logging.getLogger('ARRAY-MP')

class Station(object):
    '''
    class to store station object from metadatafile
    '''
    def __init__(self, net, sta, loc, comp,lat=0,lon=0,ele=0,dip=0,azi=0,gain=0,takeoff=0,backazi=0,sazi=0):
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
        return self.net+'.'+self.sta+'.'+self.loc+'.'+self.comp

    def pyr_name(self):
        return self.net+'.'+self.sta+'.'+self.loc+'.'

    def getcmpName(self):
        if self.loc == '--':
            self.loc=''
        return self.net+'.'+self.sta+'.'+self.loc+'.'+self.comp

    def __str__(self):
        return ('%s.%s.%s.%s')%(self.net,self.sta,self.loc,self.comp)

    def __eq__(self, other):
        return self.getName() == other.getName()

class Event(object):
    '''
    class to store event object form origin file
    '''
    def __init__(self,lat,lon,depth,time='',region='',strike=0,dip=0,rake=0):
        self.lat = lat
        self.lon = lon
        self.depth = depth
        self.region = region
        self.time = time
        self.strike= strike
        self.dip = dip
        self.rake = rake

class Trigger(object):
    '''
    class to store triggertimes for station used in xcorrelations
    '''
    def __init__(self,stationname,triggertime,arrayname,tdiff=0):
        self.sname = stationname
        self.ttime = triggertime
        self.aname = arrayname
        self.tdiff = tdiff


class Config(object):
    '''
    class to parse origin, config and metadata file for arraytool
    initialize with eventpath
    '''
    def __init__(self,eventpath):
        self.eventpath = eventpath


    def parseConfig(self,suffix):
        '''
        method to parse config files (origin,config)
        return Configdictionary
        '''
        cDict = {}
        files  = glob.glob(os.path.join(self.eventpath,'*.'+suffix))
        parser = SafeConfigParser()

        parser.read(files[0])

        for section_name in parser.sections():
            for name, value in parser.items(section_name):
                cDict[name]=value

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
                    evfile = os.path.join(self.eventpath,i)
                    fobj = open(evfile,'r')
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
                        if fnmatch.fnmatch(comp,'*HZ'):
                            MetaL.append(Station(net,sta,loc,comp,lat,lon,ele,dip,azi,gain))

            logger.info('\033[31m %d ENTRIES IN METAFILE FOUND \033[0m \n' % (len(MetaL)))
        except:
            logger.info('\033[31m METAFILE NOT READABLE \033[0m \n')

        FML = self.checkMetaInfoFile(MetaL)

        return FML

    def readpyrockostations(self):
        stations = model.load_stations(self.eventpath+'/data/stations.txt')
        MetaL = []
        for sl in stations:
                channel = sl.channels[0]
                MetaL.append(Station(str(sl.network),str(sl.station),
                str(sl.location),str(sl.channels[0])[:3],str(sl.lat),str(sl.lon),
                str(sl.elevation),str(channel.dip),str(channel.azimuth),
                str(channel.gain)))

        FML = self.checkMetaInfoFile(MetaL)

        return FML

    def readcolosseostations(self, scenario_path):
        stations = model.load_stations(scenario_path+'/meta/stations.txt')

        MetaL = []
        for sl in stations:
                channel = sl.channels[2]
                MetaL.append(Station(str(sl.network),str(sl.station),
                str(sl.location),str(sl.channels[2])[:3],str(sl.lat),str(sl.lon),
                str(sl.elevation),str(channel.dip),str(channel.azimuth),
                str(channel.gain)))
        FML = self.checkMetaInfoFile(MetaL)

        return FML




    def checkMetaInfoFile(self,MetaList):

        ML = []
        DL = []
        LL = []

        for i in MetaList:
            try:
            	if float(i.gain) == 0:
                	print 'GAIN IS ZERO ',i
                	search = ('%s.%s.%s')%(i.net,i.sta,i.loc)
               	 	DL.append(search)
               	 	LL.append(i)
    	    except:
		i.gain = (np.float(i.gain[:-3])) #careful, there is something off with some japanese/chinese stats.
                search = ('%s.%s.%s')%(i.net,i.sta,i.loc)
                DL.append(search)
                LL.append(i)


        if len(DL) > 0:
            for i in DL:
                for j in MetaList:
                    metaname = ('%s.%s.%s')%(j.net,j.sta,j.loc)
                    if i != metaname:
                        ML.append(j)
        else:
            ML = MetaList

        print len(MetaList)
        print len(ML)

        return ML


    def createFolder(self):
        '''
        method to create work folder in event directory
        returns Folderdictionary
        '''
        Folder = {}
        logger.info('\033[31m Create working environment \033[0m \n')
        if os.access(os.getcwd(), os.W_OK):

            basedir = os.path.join(self.eventpath,'work')
            sembdir = os.path.join(basedir, 'semblance')
            ascdir  = os.path.join(basedir, 'asc')
            mseeddir  = os.path.join(basedir, 'mseed')

            Folder['base'] = basedir
            Folder['semb'] = sembdir
            Folder['asc']  = ascdir
            Folder['mseed'] = mseeddir

            for key in Folder:
                if os.access(Folder[key],os.F_OK) == False:
                    os.makedirs(Folder[key])

        else:
            print "no write permissions"

        Folder['config'] = os.path.join('..','skeleton')
        return Folder


    def cpSkeleton(self,FolderDict,ConfigDict):
        '''
        method to copy skeleton scripts (plotting) to event folder
        '''
        logger.info('\033[31m Copy plotting scripts \033[0m \n')

        cpf = ConfigDict['plotcopy'].split(',')
        for i in cpf:
            if os.access(os.getcwd(), os.W_OK):
                try:
                    shutil.copy2(os.path.join(FolderDict['config'],ConfigDict[i]),os.path.join(FolderDict['semb'],ConfigDict[i]))
                except:
                    logger.info('\033[31m Copy processing scripts Error \033[0m \n')
                    continue


    def writeConfig(self,Config,Origin,Folder):
        '''
        method to write recently used config to event folder
        '''
        fobj = open(os.path.join(Folder['semb'],'stations_0.cfg'),'w')
        fobj.write('%source: '+Origin['lat']+' '+Origin['lon']+' '+Origin['depth']+' '+Origin['time']+'\n')
        for i in Config:
            fobj.write('%'+i+': '+Config[i]+'\n')
        fobj.close()


    def writeStationFile(self,StationMetaData,Folder,flag):
        '''
        method to write recently used stations for processing to event folder
        '''
        name = 'stations_'+str(flag)+'.dat'
        fobj = open(os.path.join(Folder['semb'],name),'w')
        for i in StationMetaData:
            fobj.write('%s %s %s %s\n' %(flag,i.getName(),i.lat,i.lon))
        fobj.close()
