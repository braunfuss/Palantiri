
'''
modul for deserializing pickle data from different processes
'''

import os
import sys
import platform

WINDOWS = (platform.system() == 'Windows')

# add local directories to import path 

sys.path.append ('../tools/')                     
sys.path.append ('../Common/')
                  
import cPickle as pickle

#       Import from common

import  Basic
import  Globals
import  Logfile
import  Debug

#       Import from Process

from ttt import MinTMaxT

# -------------------------------------------------------------------------------------------------

def deserializeTTT (numproc):

        L = []

        for i in range(numproc) :
            data = Basic.loadDump (str(i)+'-ttt.pkl')

            if data != None : L.append(data)
        #endfor
    
        TTTGridMap = {}

        for i in L:
            for j in i.iterkeys() : TTTGridMap[j] = i[j]

        return TTTGridMap

# -------------------------------------------------------------------------------------------------

def deserializeMinTMaxT (numproc):

        L = []

        for i in range(numproc):
            data = Basic.loadDump ('minmax-'+str(i)+'.pkl')

            if data != None : L.append(data)
        #endfor

        mint = min ([x.mint for x in L])
        maxt = max ([x.maxt for x in L])

#        for i in range(numproc):
 #           os.remove ('minmax-'+str(i)+'.pkl')
        
        return mint,maxt

# -------------------------------------------------------------------------------------------------

def deserializeSembDict(numproc):

        L = []

        for i in range(numproc):
            data = Basic.loadDump ('sembDict-'+str(i)+'.pkl')

            if data != None : L.append(data)
        #endfor

        sembDict = {}

        for i in L:
            for j in i.iterkeys() : sembDict[j] = i[j]
            
        return sembDict
    
# -------------------------------------------------------------------------------------------------

def deserializeSemb_unused (numproc):
            print 'ENTER DES SEMB\n'

            L = []
        #for i in range(numproc):
            fname = '1-semb.pkl'
            print fname
            data = Basic.loadDump (fname)

            L.append(data)

            sembDict = {}
            print L

        #for i in L:
         #   for j in i.iterkeys():
          #      sembDict[j] = i[j]
            
    #    for i in range(numproc):
     #       os.remove('sembDict-'+str(i)+'.pkl')

            return L
# -------------------------------------------------------------------------------------------------
    
def deserializeSembMaxFile (numproc):

    L = []

    for i in range(numproc):
        data = Basic.loadDump ('sembMAX-'+str(i)+'.pkl')

        if data != None :  L.append(data)
    #endfor

    sembMax = {}

    for i in L:
        for j in i.iterkeys() : sembMax[j] = i[j]
            
    return sembMax
