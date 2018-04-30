import sys

from   Program  import startTest           # from Common
import Basic

Basic.onlyForInternalUse ()

sys.argv = startTest ('process', workingDir='tmpProcess')

import main                                # from Process
main.MainProc ()

 