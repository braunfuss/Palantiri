import os
import fnmatch
import sys

def findOriginFile(event):
    
    evfile = ''
    path = os.path.join(os.getcwd(),'..','events',sys.argv[2])
    for i in os.listdir(path):
        if fnmatch.fnmatch(i, event+'*.origin'):
                evfile = os.path.join(path,i)

    return evfile