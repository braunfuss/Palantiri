
import os
import sys

#   add local directories to import path
                   
sys.path.append('Common/')

#   import from Common

#import Basic

# -------------------------------------------------------------------------------------------------

def start() :

    if sys.argv[1] != 'new_version' : return False
       
    #  python arraytool.py(0)  new_version(1) 
         
    at = os.path.join(os.getcwd(),'Common', 'NewVersion.py')      # python script

    os.system(sys.executable + ' ' + at)
    return True


