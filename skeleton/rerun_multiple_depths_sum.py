import os
import sys
import fnmatch
import shutil
from optparse import OptionParser
import logging

logger = logging.getLogger(sys.argv[0])
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

def getDepths():
    L=[]
    for i in os.listdir(os.getcwd()):
       if fnmatch.fnmatch(i,'*_*.ASC'):
            i = i.split('_')
            L.append(i[0])
    K = list(set(L))

    K= sorted(K)
    return K

def run(begin,end,depths):

    fd = open('duration.txt','r')
    fw = open('sembmaxvalue.txt','r')

    fdc=0
    for i in fd:
        if fdc == 0:
            line1 = i
        fdc+=1    
    fd.close()

    c=0
    for i in fw:
        line = str.split(i,' ')
        if fnmatch.fnmatch(line[1], begin):
            beginline = c
            begincontext = line
        if fnmatch.fnmatch(line[1], end):
            endline = c
            endcontext = line
        c+=1    

    fw.seek(0)
        
    fnd = open('tmp.txt','w')
    fnd.write(line1)
    line2 ='# Shifted event onset at i = '+ begincontext[0] +' time ='+begincontext[1]+' semblance = '+begincontext[2]+'\n'
    fnd.write(line2)

    counter=0
    for i in fw:
            if counter >= beginline and counter <=endline:
                fnd.write(i)
                line = str.split(i)
            counter+=1
    
    line3='# Shiftet event stop at i = '+endcontext[0]+' time = '+endcontext[1]+' semblance = '+endcontext[2]+'\n'
    fnd.write(line3)
    duration = int(end)-int(begin)
    ll ='# Event duration = '+str(duration)
    fnd.write(ll)
    fnd.close()

    shutil.copyfile('tmp.txt', 'duration.txt')

    D=[]
    fobj = open('duration.txt','r')
    for line in fobj:
            if line[0] != '#':
                line = line.split()
                D.append(float(line[2]))
    fobj.close()
    minVal= min(D)
    maxVal= max(D)

    fobj = open('minmax.txt','w')
    fobj.write(('%f,%f')%(minVal,maxVal))
    fobj.close()

    cmd='./summary-plot-point_multiple_depths.sh test '+depths
    print cmd
    os.system(cmd)

def main(args):
    parser = OptionParser(usage="%prog -s 175 -e 190 -d 15.5")
    parser.add_option("-s", "--start", type="string", dest="begin", help="begin in s")
    parser.add_option("-e", "--end", type="string", dest="end", help="end in s")
    parser.add_option("-d", "--depth", type="string", dest="depth", help="depth in km")

    return parser.parse_args(args)

if __name__ == "__main__":
    D = getDepths()
    
    options,args = main(sys.argv)
    
    print options.begin,options.end,options.depth,type(options.begin),type(options.end),type(options.depth)
    
    if len(sys.argv) == 1: 
        os.system(('%s %s --help')%(sys.executable,sys.argv[0]))
        sys.exit()
    
    if not options.begin:
        os.system(('%s %s --help')%(sys.executable,sys.argv[0]))
        logger.info('\033[31mStart value duration is missing \033[0m')
        sys.exit()
    if not options.end:
        os.system(('%s %s --help')%(sys.executable,sys.argv[0]))
        logger.info('\033[31mEnd value duration is missing \033[0m')
        sys.exit()
    if not options.depth:
        os.system(('%s %s --help')%(sys.executable,sys.argv[0]))
        logger.info('\033[31mDepth value is missing \033[0m')
        logger.info('\033[31mAvailable Depths: \033[0m')
        for i in D:
            logger.info('\033[31m%s km\033[0m'%(i))
        sys.exit()
    
    
    run(options.begin,options.end,options.depth)

