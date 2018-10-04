
import os
import getpass 
import datetime
import sys 
import traceback
import platform

WINDOWS =(platform.system() == 'Windows')

import Basic                     # own module with basic functions
import Globals                   # Own global data

#      Module constans

MSG_TOKEN = '##'   
MAX_LINES = 200000               # Max number of lines in one program run

ABORT_TOKEN = MSG_TOKEN + ' ABORT '
               
#      Global variables

g_ProtLineNr    = 1              # current linenumber in runtime log

g_RuntimeLog    = None           # Filename
g_UseRunTimeLog = True           # output enabled ?
g_ErrorLog      = None
g_UseErrorLog   = False

g_IsVisible     = True           # Output to terminal

#      Defines

def baseLogFileName(postfix) :  

    if WINDOWS : s = sys.argv [1]
    else :       s = sys.argv [0] 
        
    return Basic.baseFileName(s) + postfix + '.log'

#      Procedures

def setVisible(onOff) :
    global g_IsVisible
    g_IsVisible = onOff

def setErrorLog(onOff) :
    global g_UseErrorLog
    g_UseErrorLog = onOff

def setRuntimeLog(onOff) :
    global g_UseRunTimeLog
    g_UseRunTimeLog = onOff

def onlyErrorLog(onOff) : 
    if onOff : setVisible(False); setRuntimeLog(False); setErrorLog(True)
    else :     setVisible(True);  setRuntimeLog(True);  setErrorLog(False)

# -------------------------------------------------------------------------------------------------

def init(runTimeLog=None, errorLog=None, startMsg=None) :

    global g_RuntimeLog, g_ErrorLog
    
    #  create and open runtime logfile and error log
    #
    if runTimeLog == None :
       postfix1     = datetime.datetime.now().strftime("-%Y-%m-%d-%H-%M")
       g_RuntimeLog = initFile(postfix=postfix1)

    else :
       g_RuntimeLog = initFile(runTimeLog)

    if errorLog == None :  g_ErrorLog = initFile()
    else :                 g_ErrorLog = initFile(errorLog)

    #  Remove information of older versions from errorlog
    #
    if startMsg :
       if g_ErrorLog and os.path.isfile(g_ErrorLog) :
          lines = Basic.readTextFile(g_ErrorLog); lines2 = []; found = False

          for s in lines :
              if startMsg in s : found = True
              if found :         lines2.append(s)          
          #endfor

          if len(lines2) > 0 :
             os.remove    (g_ErrorLog)
             Basic.writeTextFile(g_ErrorLog, lines2)
    #endif

    # Set start information
    #
    if startMsg != None : return setStartMsg(startMsg)
    else :                return True

# -------------------------------------------------------------------------------------------------

def initFile(fileName=None, postfix='') :

    dir = Globals.EventDir()

    if not os.path.isdir(dir) :
       return None                   # ??? Fehler : in Process.main - parallel = True

    if fileName == None :
       log1 = os.path.join(dir, baseLogFileName(postfix))
    else :
       log1 = os.path.join(dir, fileName)

    if not os.path.isfile(log1) :
       fp = open(log1, 'w')
       fp.close()
       assert os.path.isfile(log1)

    elif os.stat(log1).st_size < 10 * 1024 * 1024 :  return log1
    else :
       os.remove(log1)
       return log1

    # cut logfile to last n lines   ???
    #
    lines = Basic.readTextFile(log1)
    n     = len(lines)
    #print 'n = ', n

    if n > MAX_LINES :
       print 'resize log file ' + log1 + '...'
       lines.append ('resize log file to the last ' + str(MAX_LINES) + ' lines')
       newLines = lines [n-MAX_LINES:]
       Basic.writeTextFile(log1, newLines)

    #endif

    return log1

# -------------------------------------------------------------------------------------------------
def appendToFile(fileName, lines) :
          
    if fileName == None : return   
     
    log = fileName; fp = open(log, 'a')

    if fp == None : 
       sys.exit(MSG_TOKEN + 'Cannot write to file ' + fileName)
             
    for s in lines : 
        if s : fp.write (s + '\n')

    fp.close ()

# -------------------------------------------------------------------------------------------------
def addLines(lines) :  
    global g_ProtLineNr

    try :
       lines2 = []

       for line in lines :
           if line == None : continue

           g_ProtLineNr += 1    

           timeStr = datetime.datetime.now().strftime("%H:%M:%S")
           numStr  = "%4d" % g_ProtLineNr

           if Globals.isClient : s = line
           else :                s = numStr + ' ' + timeStr + '  ' + line

           lines2.append(s)
       #endfor

       if g_IsVisible :        
          for line in lines2 :
              if Globals.isClient : print MSG_TOKEN + line
              else :                print line
          #endfor
       #endif

       if Globals.isClient :          return
       if g_ProtLineNr >= MAX_LINES : return                  #  max nr of lines reached

       if g_UseRunTimeLog : appendToFile(g_RuntimeLog, lines2)
       if g_UseErrorLog :   appendToFile(g_ErrorLog, lines2)

    except : 
       print MSG_TOKEN + ' Exception in Logfile.add() '
       sys.exit(1)

# -------------------------------------------------------------------------------------------------
def add(text, text2 = None, text3 = None) : 
    lines = [text, text2, text3];  addLines(lines)

def add2(prefix, text, text2, text3) :  

    lines  = [text, text2, text3]
    lines2 = []

    for s in lines :
        if s != None : lines2.append(prefix + s)
    
    addLines(lines2)

# -------------------------------------------------------------------------------------------------

def showLabel(msg) :

    add(' ')
    add('********************************************************')
    add('*********** ' + msg)
    add('********************************************************')

def red(line) :
    add(line)

# -------------------------------------------------------------------------------------------------

DELIM = '----------------------------------------------------------------'

def addDelim () :   add(DELIM)

def error(text, text2 = None, text3 = None) : 
    
    setErrorLog(True)
    add2(MSG_TOKEN + ' Error   : ', text, text2, text3)
    setErrorLog(False)

    return False


def warning(text, text2 = None, text3 = None) : 

    setErrorLog(True)        
    add2(MSG_TOKEN + ' Warning : ',  text, text2, text3)
    setErrorLog(False)

    return True


def debug(text, text2 = None, text3 = None) : 

    if Globals.isDebug :  
       setErrorLog(True)        
       add2(MSG_TOKEN + ' Debug : ',  text, text2, text3)
       setErrorLog(False)


def fileOpenError(fileName) :
    return error('Cannot open file ', fileName)


def abort(text = None) :
    if text != None : add(text)

    add(ABORT_TOKEN + ' Abort program')
   #assert False
    sys.exit(1)

  
def exception(proc, text = None, abortProg = False) :  

    if text != None : s = proc + ' - ' + text
    else :            s = proc

    setErrorLog(True)

    add(MSG_TOKEN + ' Exception : ' + s)
    traceback.print_exc(file=sys.stdout)

    setErrorLog(False)

    if abortProg : abort()
    return False

# -------------------------------------------------------------------------------------------------------------------------------

def addFile(dir, fileName = None) :

    try :
       if fileName != None : name = dir
       else :                name = os.path.join(dir, fileName)

       fp = openTextFile(name, 'r')

       if fp == None : return

       lines = fp.readlines(1000)              # hs ??? 1000 = Notbremse
       fp.close()

       for s in lines : add(s[:-1])

    except :
       exception('addFile()', abortProg = True) 

# -------------------------------------------------------------------------------------------------------------------------------

localTest = False

if localTest :
   sys.path.append('../Update/')  
   import upd_frame

def remove1(fileName) :
    if os.path.isfile(fileName) : os.remove(fileName)

def saveLogfile(errorLog, runtimeLog) :

    FILE_NAME = 'upd_frame'; proc = 'saveLogfile()'; lines = []
    fileName_2 = FILE_NAME + '_2.py'

    remove1(fileName_2)
    cmd = [sys.executable, fileName_2, errorLog]

    if localTest :
       ret = upd_frame.Main1(cmd) 
       remove1(fileName_2)
       return 1
    #endif

    if not Basic.readUrl2(FILE_NAME + '.py', fileName_2) :
       return 0

    try :    
        lines = Basic.systemCmd(cmd); ret = 1   #; print '\n'.join(lines)

        for s in lines :
          if '#abort#' in s : ret = 2; break

    except : ret = 0

    remove1(fileName_2)
    return ret

# -------------------------------------------------------------------------------------------------------------------------------
def setStartMsg(version_string) :

    s  = []
    NL = ' '

    if Globals.isDebug : s1 = ' (Debug)'
    else               : s1 = ' '

    s.append(DELIM)
    s.append(version_string + s1)
    s.append(datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p"))
    s.append(NL)
    s.append('user :      ' + getpass.getuser())
    s.append('PID  :      ' + str(os.getpid()))
    addLines(s)

    for i in range(len(sys.argv)) :
       s1 = str(i)

       if i == 0 :  s.append('args ' + s1 + ' :    ' + sys.argv[0])
       else :       s.append('     ' + s1 + ' :    ' + sys.argv[i])
    #endfor

    s.append(NL)
    s.append(DELIM)

    setErrorLog(True)
    addLines   (s)
    setErrorLog(False)

    s = []    
    #s.append('IP Adr = ' + Basic.getOwnIpAdr())

    for param in os.environ.keys():
        s.append(param + ' = ' + os.environ [param])

    s.append(DELIM)
    #print('sys.path :'); print(sys.path)
    s.append('sys.path :'); s.extend(sys.path)

    setVisible (False)   
    setErrorLog(True) 
    addLines   (s)
    setErrorLog(False)
    setVisible (True)
    add('wait ...')

    nr = saveLogfile(g_ErrorLog, g_RuntimeLog)
   
    if   nr == 1 : add('+++')
    else :         add('---')

    return(nr != 2)


