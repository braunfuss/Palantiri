
#      Import from tools

import Basic                     # own module with basic functions
import Globals                   # Own global data
import Logfile

#      Module constans

#      Global variables

g_Level = 0
g_showProcs = True

# -------------------------------------------------------------------------------------------------
#
#     class PROC implements enter/leave protocol              ??? funktioniert noch nicht
#
class PROC(object):

   def __init__(self, name):
      global g_Level
      assert g_showProcs != None

      self.name = name
      self._print1('E')
      g_Level += 1

   def __del__(self):
      global g_Level

      g_Level -= 1
      self._print1('L')

   def _print1(self, text):
       global g_showProcs

       xxx
       if not g_showProcs: return

       s = ''

       for i in range(g_Level): s.append('  ')


#endClass PROC

# -------------------------------------------------------------------------------------------------

def assert1(condition, text1 = None, text2 = None, text3 = None):

    if not condition:
       if text1 != None: Logfile.error(text1,text2,text3)

       Logfile.abort()
    #endif

def assert2(condition1, condition2,  text1 = None, text2 = None, text3 = None):

    assert1(condition1, text1,text2,text3)
    assert1(condition2, text1,text2,text3)

# -------------------------------------------------------------------------------------------------

def enter(module, procName):
    global  g_Level

    if not Globals.isDebug: return

    s = ''
    for i in range(g_Level): s.append('  ')

    Logfile.add(s + 'E ' + module + '.' + procName)
    g_Level += 1


def leave(module, procName, val=0):
    global  g_Level

    if not Globals.isDebug: return

    g_Level -= 1

    s = ''
    for i in range(g_Level): s.append('  ')

    Logfile.add('L ' + module + '.' + procName)


def init():
    global  g_showProcs

    g_showProcs = True
    return

    key = 'debug_enter'   # ???

    if key in os.environ.keys(): g_showProcs =(int(os.environ [key]) > 0)
    else:                        g_showProcs = False
