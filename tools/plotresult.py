import logging
from optparse import OptionParser
import os,sys
import subprocess
import fnmatch
import config

logger = logging.getLogger('ARRAY-MP')
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

def init():
    parser = OptionParser(usage="%prog -f Eventpath ")
    parser.add_option("-f", "--evpath", type="string", dest="evpath", help="evpath")

    (options, args) = parser.parse_args()
    if options.evpath == None:
        parser.error("non existing eventpath")
    evpath = options.evpath
    plotpath = os.path.join(evpath,'work','semblance')
    
    return plotpath,evpath

def plotresults(pathtoplot,eventpath):


    C = config.Config(eventpath)
    Config = C.parseConfig('config')
    ntimes = int((int(Config['forerun']) + int(Config['duration']))/int(Config['step']))
    
    os.chdir(pathtoplot)
    
    cmd = ('%s plot.py')%(sys.executable)
    
    L = []
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for index,line in enumerate(iter(p.stdout.readline, "")):
        L.append(line)
        if fnmatch.fnmatch(line,'*_semblance_*'):
            sys.stdout.write(('\033[31m\rArrayTool plot progress: timestep  %d of %s\033[0m')%(index/3,ntimes))
        if fnmatch.fnmatch(line,'summary_*'):
            sys.stdout.write(('\033[31m\rArrayTool summary plot: %s\033[0m')%(line))
        sys.stdout.flush()
    p.wait()
    

    sys.stdout.write(('\033[31m\rArrayTool plot results finished\033[0m')%())
    sys.stdout.write(('\033[31m\rArrayTool plot results saved in %s \n\033[0m')%(pathtoplot))


if __name__ == "__main__":
    options,epath = init()
    
    plotresults(options,epath)