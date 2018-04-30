
import os
import sys
import Basic

#VERSION_STRING = 'Version 2.0 - 18.September 2014'
VERSION_STRING  = 'Version 2.0 - 14.November 2014'

def versionLine() : return 'VERSION_STRING = ' + '\'' + VERSION_STRING + '\''

DIRECTORIES = ['Cluster','Process','Waveform']

# -------------------------------------------------------------------------------------------------
# Generate all new files 'Version.py' for all 'directories'

def generate () :

    user = Basic.getUserName()

    if user != 'staedtke' : return

    print '\n' + VERSION_STRING + '\n' + user + '\n'

    for dir in DIRECTORIES :
        file = os.path.join (os.getcwd(),dir, 'Version.py')
        print 'file = ', file

        fp = open (file, 'w')

        if fp == None :
           print 'Cannot open file ' + file
           continue

        fp.write (versionLine() + '\n')
        fp.close ()
    #endfor

    print ' '
    print '********************************************************'
    print '*********** New version generated '
    print '********************************************************\n'

# -------------------------------------------------------------------------------------------------
# Check if all files 'Version.py' have the same contents

def check() :

    lines = []

    for dir in DIRECTORIES :
        file = os.path.join (os.getcwd(),dir, 'Version.py')

        fp = open (file, 'r')

        if fp == None : 
           print 'Cannot open file ' + file
           sys.exit (1)
   
        line = fp.readline()
        lines.append (line[:-1])
        fp.close ()
    #endfor

# -------------------------------------------------------------------------------------------------    
        
if __name__ == "__main__" :

    generate ()

    

