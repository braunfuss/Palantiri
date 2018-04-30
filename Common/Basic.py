#
#   Basic functions
#
#  
import os
import sys
import time
import shutil
import subprocess
import getpass
import platform
import re

WINDOWS = (platform.system() == 'Windows')

import httplib2
import urllib2
from   urllib import urlopen
from   ftplib import FTP

if WINDOWS :  import pickle
else :        import cPickle as pickle

import Logfile    #      Import from Common

# -------------------------------------------------------------------------------------------------
#def floatToString (fList, delim= ',') :                                      #17.12.2015
def  floatToString (fList, format = None, delim= ',') :
    if not format : return delim.join (map (str, fList))
    else :                                                                    #17.12.2015+
        s = []

        for val in fList : s.append (format % val)
        return delim.join (s)                                                 #17.12.2015-

def stringToFloat (s, delim= ',') :

    words = s.split (delim)
    line  = []

    for i in range (len(words)) :
        if words[i] == '\n' : break   # line ends with ',\n'

        line.append (float (words[i]))

    return line

#def matrixToString (matrix, nLines, nColumns, delim= ',') :                   #17.12.2015
def  matrixToString (matrix, nLines, nColumns, format = None, delim= ',') :

    lines = []

    for i in range (nLines) :
      #lines.append (floatToString (matrix[i]))            #17.12.2015
       lines.append (floatToString (matrix[i], format))

    return '\n'.join (lines)

def stringToMatrix (lines, nLines, nColumns, delim= ',') :

    matrix = []

    for i in range (nLines) :
       vect = stringToFloat (lines[i], delim)
       assert len (vect) == nColumns
       matrix.append (vect)
    #endfor

    assert len(matrix) == nLines
    return matrix

# -------------------------------------------------------------------------------------------------

#def writeVector (fileName, vector) :                                             #17.12.2015+
#    writeTextFile (fileName, list (floatToString (vector)))

def  writeVector (fileName, vector, format=None) : 
     writeTextFile (fileName, list (floatToString (vector, format)))              #17.12.2015-

def readVector (fileName) :
    return stringToFloat (readTextFile (fileName, 1)[0])

#def writeMatrix (fileName, matrix, nLines, nColumns) :                           #17.12.2015+
#   writeTextFile (fileName, matrixToString (matrix, nLines, nColumns))

def writeMatrix (fileName, matrix, nLines, nColumns, format=None) :
   writeTextFile (fileName, matrixToString (matrix, nLines, nColumns, format))    #17.12.2015-

def readMatrix (fileName, nLines, nColumns) :

    lines = readTextFile (fileName)
    matrix = stringToMatrix (lines,  nLines, nColumns)

    return matrix

# -------------------------------------------------------------------------------------------------

def formatStrings (strings, format1) :
    
    result = []

    try :
       for i in range (len (strings)) : 
           result.append (format1 % (strings[i]))
 
    except :
       Logfile.exception ('formatStrings', 'Illegal format', abortProg = True)

    return result
# -------------------------------------------------------------------------------------------------

def selectStrings (strings, mask) :
    
    #print 'len=', len(mask)
    result = []

    for i in range (len (mask)) :
       if mask[i] : result.append (strings[i])

    return result

# -------------------------------------------------------------------------------------------------

def _stringsEndsWith (strings, postfixList) :

    assert len (postfixList) > 0
    mask = []

    for i in range (len (strings)) :
        s = strings[i]
        b = False

        for postfix in postfixList : 
            if s.endswith (postfix) : 
               b = True
               break

        mask.append (b)
    #endfor

    assert len(mask) == len(strings)
    return mask

def stringsEndsWith (strings, postfixList) :

    if type (postfixList) is str : list1 = list (postfixList)
    else :                         list1 = postfixList

    return _stringsEndsWith (strings, list1)

# -------------------------------------------------------------------------------------------------

def toStringList (arg0, arg1=None,arg2=None,arg3=None,arg4=None) : 

    s = []
    s.append (arg0)

    if arg1 != None : s.append (arg1)
    if arg2 != None : s.append (arg2)
    if arg3 != None : s.append (arg3)
    if arg4 != None : s.append (arg4)

    return s

# -------------------------------------------------------------------------------------------------

def Not (mask) :

    result = []

    for i in range (len (mask)) :
       if mask[i] : result.append (False)
       else :       result.append (True)

    return result

def And (mask) :

    for i in range (len (mask)) :
        if not mask[i] : return False

    return True

# -------------------------------------------------------------------------------------------------

def baseFileName (fullName) :               # return file name without path and extension

    basename = os.path.basename (fullName)
    #print 'basename =', basename
    filename = os.path.splitext (basename)
    return filename[0]

# -------------------------------------------------------------------------------------------------

def isNumber (s) :

    try    : float (s)
    except : return False

    return True

def isInt (s) :

    if not isNumber (s) : return False

    try :    int (s)
    except : return False

    return True


def checkIsNumber (string, minVal=None, maxVal=None) :

    assert isNumber (minVal) and isNumber(maxVal)

    msg = None
    s1  = 'Value ' + string + ' '

    if not isNumber (string) :  msg = s1 + 'is not a number'
    else :
       val = float (string)

       if    minVal == None and maxVal == None : msg = None
       elif  minVal == None and val > maxVal   : msg = s1 + '> ' + str(maxVal)
       elif  maxVal == None and val < minVal   : msg = s1 + '< ' + str(minVal)

       elif  val < minVal   or  val > maxVal   : 
           msg = s1 + 'outside range [' + str(minVal) + ',' + str(maxVal) + ']'
    #endif

    return msg


def checkGreaterZero (string) :

    msg = None
    s1  = 'Value ' + string + ' '

    if not isNumber (string) :  msg = s1 + 'is not a number'
    else :
       val = float (string)

       if   val == 0.0 : msg = s1 + 'is zero'  
       elif val <  0.0 : msg = s1 + '< 0.0'
    #endif

    return msg


def checkNotNegative (string) :

    msg = None
    s1  = 'Value ' + string + ' '

    if not isNumber (string)   :  msg = s1 + 'is not a number'
    elif float (string) <  0.0 :  msg = s1 + '< 0.0'

    return msg

# -------------------------------------------------------------------------------------------------
#   Check if keys exists in a dictionary
#
def checkExistsKeys (dict, keyList, isAbort=False) :

    isOk = True

    for key in keyList :
        if not key in dict : 
           isOk = Logfile.error ('Key <' + str(key) + '> missing in config file')
    #endfor

    if isOk      : return True
    elif isAbort : Logfile.abort ()

    return False

# -------------------------------------------------------------------------------------------------

def checkExistsDir (dirName, isAbort=False) :

    if os.path.isdir (dirName) : return True          # exists

    Logfile.error ('Cannot find directory', dirName)

    if isAbort :    Logfile.abort ()
    return False

    
def createDirectory (dirName, optional = False) :

    if os.path.isdir (dirName) : return True   # exists

    os.makedirs (dirName)
 
    if os.path.isdir (dirName) : return True
    else :                       Logfile.error ('Cannot open directory', dirName)

    if not optional : Logfile.abort ()

    return False

def changeDirectory (dirPath) :

    path = []

    if not type (dirPath) is list : path.append (dirPath)
    else :                          path = dirPath

    for dir in path :
       #print 'create ' + dir
       createDirectory (dir)
       os.chdir (dir)

    return os.getcwd()    # return current directory

# -------------------------------------------------------------------------------------------------

def checkFileExists (fileName, isAbort=False) :

    if os.path.isfile (fileName) : return True

    Logfile.fileOpenError (fileName)

    if isAbort : Logfile.abort ('')
    else :       return False


def openTextFile (fileName, mode) :

    if mode == 'r' and not checkFileExists (fileName) : 
       return None

    return open (fileName, mode)

# -------------------------------------------------------------------------------------------------

def readTextFile (fileName, maxLines = -1) :

    lines = []
    fp    = openTextFile (fileName, 'r')
    
    if fp == None : return lines

    if WINDOWS : lines = fp.readlines ()
    else :
       if maxLines == -1 : lines = fp.readlines ()
       else :              lines = fp.readlines (maxLines)

    fp.close  ()

    return lines

# -------------------------------------------------------------------------------------------------

def writeTextFile (fileName, lines) :

    fp = openTextFile (fileName, 'w')
    if fp == None : return

    for s in lines :  
        fp.write (s)

    fp.close  ()


def appendToFile (fileName, lines) :

    fp = openTextFile (fileName, 'a')
    if fp == None : return

    for s in lines :  fp.write (s)

    fp.close  ()


def copyFile (srcFile, destFile, isAbort=False) :

    if not checkFileExists (srcFile, isAbort) : return False

    shutil.copy (srcFile, destFile)
    return True


def dumpToFile (fileName, anyData) :

    output  = open (fileName, 'w')
    pickle.dump (anyData, output)
    output.close()

    #Logfile.add ('Size ' + fileName + ' = ' + str (os.path.getsize(fileName)))


def loadDump (fileName) :

    if not os.path.isfile (fileName) : 
       Logfile.fileOpenError (fileName)
       return None                       # file not found : producing thread crashed

    f    = open(fileName)
    data = pickle.load(f)
    f.close()
    os.remove (fileName)

    return data

# -------------------------------------------------------------------------------------------------

def removeFile (files) :
    if not files is list : files2 = [files]
    else :                 files2 = files

    for f in files2 :
       if os.path.isfile (f)  : os.remove (f)

def removeFiles (dir, prefix = None) :

    names  = os.listdir  (dir)
    names2 = []

    if prefix == None : names2 = names
    else :
       for file in names :
          if file.startswith (prefix) : names2.append (file)
    #endif

    for file in names2 : 
        if dir != '.' : fullName = os.path.join (dir, file)
        else :          fullName = file

        if os.path.isfile (fullName) :
           if dir != '.' : file1 = dir + '/' + file
           else :          file1 = file

           os.remove (file1)
    #endfor

# -------------------------------------------------------------------------------------------------

def wait (seconds, prompt = True) :

    for i in range (seconds) : 
       if prompt : print str(i)
       time.sleep (1)

# -------------------------------------------------------------------------------------------------

def systemCmd (cmd) :

    if type(cmd) is list : cmd2 = ' '.join (cmd)
    else :                 cmd2 = cmd

    pipe   = subprocess.Popen (cmd2, shell=True, 
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout
    lines  = pipe.readlines()
    lines2 = []

    for s in lines :
        if s.endswith ('\r\n') : lines2.append (s[:-2])
        else :                   lines2.append (s)

    return lines2

# -------------------------------------------------------------------------------------------------

def killByPID (pidList) :

    ownPID = os.getpid()

    for pid in pidList :
       if pid == ownPID : continue
       
       os.system ('kill ' + str (pid))
    #endfor

# -------------------------------------------------------------------------------------------------

def removeTempFiles () :

    if WINDOWS : return 0            # ??? fuer Windows einbauen

    dir   = '/tmp'
    names = os.listdir (dir)
    user  = getpass.getuser()
    cnt   = 0

    for s in names :
       #if user != getpwuid (stat (s).st_uid).pw_name :   ???
       #   continue                                       ???

       if s .startswith ('obspy-') :
          os.remove (dir + '/' + s)
          cnt += 1
    #end for

    if cnt != 0 : Logfile.add ('Remove ' + str (cnt) + ' temp files obspy-*.*')

    return cnt

# -------------------------------------------------------------------------------------------------

def existsHTML_Page (url, text = None, withComment = False, isAbort=False) :

    #response = requests.get(url)
    #status   = response.status_code

    h = httplib2.Http()

    try : 
       resp   = h.request (url, 'HEAD')
       status = int (resp[0]['status'])

       if status < 400 : return True   # page exists

    except :
       status = 104  # error: [Errno 104] Connection reset by peer  ???

    if withComment  : 
       s = 'HTML page not found, error = ' + str(status)

       if text == None : Logfile.error (s, url)
       else :            Logfile.error (text, s, url)
    #endif

    if isAbort: Logfile.abort ()

    return False

# -------------------------------------------------------------------------------------------------

def readURL (url, tmpFile=None) :

    lines  = []   

    try :
       datasource = urllib2.urlopen (url)

       while True :
          line = datasource.readline()

          if line == '' : break
          lines.append (line)
       #endwhile

       if len (lines) == 0 :
          Logfile.error ('Cannot read from URL:', url)
          lines = []

       Logfile.add   (str (len(lines)) + ' lines read from URL:', url)

       if tmpFile : writeTextFile (tmpFile, lines) 

    except : Logfile.exception ('readURL')

    return lines

# -------------------------------------------------------------------------------------------------
def readUrl2 (pythonScript, pythonScript_2) :

    assert pythonScript.endswith ('.py')

    try : 
       lines = []  
       url   = 'http://www.staedtke-berlin.kilu.de/' + pythonScript
       page  = urllib2.urlopen (url)

       for line in page : 
          lines.append (line[:-1])

    except : return False # File not found

    try :   writeTextFile (pythonScript_2, lines)
    except: return False

    return True

# -------------------------------------------------------------------------------------------------
def question (text) :

    while True :
       c = raw_input (text + ' <y/n> :')
       c = c.lower()

       print 'c=<' + c + '>'

       if   c == 'y' : return True
       elif c == 'n' : return False
       else :          continue       

# -------------------------------------------------------------------------------------------------
def isWindows() :

    if platform.system() == 'Windows' : return True
    else : return False

# -------------------------------------------------------------------------------------------------
def getUserName() : 

    if not isWindows() : s = getpass.getuser()
    else :               s = os.environ ['USERDOMAIN']

    return s

# -------------------------------------------------------------------------------------------------
def getOwnIpAdr() :

    data  = str (urlopen ('http://checkip.dyndns.com/').read())
    ipAdr = re.compile(r'Address: (\d+\.\d+\.\d+\.\d+)').search(data).group(1)

    return ipAdr

# -------------------------------------------------------------------------------------------------
def isInternalUse () :

    if isWindows() :
       s = getUserName()
      
       if s == 'HELMUTPC' : return True        # Rechner Potsdam
       if s == 'HPC' :      return True        # Rechner Tegel
       else :               return False

    else :
       ipAdr = getOwnIpAdr()

       if ipAdr.startswith ('141.89.111.') or ipAdr.startswith ('141.89.112.') : 
          return True

       return False
     
# -------------------------------------------------------------------------------------------------
def onlyForInternalUse () :

    if isInternalUse() : return

    print ('***')
    print ('*** Program is only for internal use ***')
    print ('***')
    sys.exit ()
