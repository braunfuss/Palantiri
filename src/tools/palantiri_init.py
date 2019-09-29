import os
import fnmatch
import sys
from optparse import OptionParser
import logging
import imp
import obspy.core
from palantiri.process import ProcessProgs
from palantiri.cluster import ClusterProgs
from palantiri.common import CommonProgs
from palantiri.tools import config as config
import palantiri
import shutil

logger = logging.getLogger(sys.argv[0])
logger.setLevel(logging.DEBUG)


def main():
    if len(sys.argv) < 2:
        print("workdir path/name missing")
        quit()

    foldername = sys.argv[1]
    workfolder = os.path.join(os.getcwd(), './', foldername)
    eventsfolder = os.path.join(os.getcwd(), './', foldername, 'events')
    tttgridsfolder = os.path.join(os.getcwd(), './', foldername, 'tttgrid')
    tmpfolder = os.path.join(os.getcwd(), './', foldername, 'tmpProcess')

    if os.access(workfolder, os.F_OK) is False:
            os.makedirs(workfolder)
            os.makedirs(tmpfolder)
            os.makedirs(tttgridsfolder)
            os.makedirs(eventsfolder)

            logger.info('\033[31m Working Super-FOLDER CREATED \033[0m \n')

    else:
        print("workdir already exists!")
        quit()


    dstfolder = foldername
    dstfile = ('global.conf')
    path = palantiri.__path__
    src = os.path.join(path[0], 'skeleton', 'global.conf')
    dst = os.path.join(dstfolder, dstfile)
    logger.info('\033[31m Created work directory \
                %s \033[0m \n' % (dstfolder.split('/')[-1]))
    shutil.copy(src, dst)
