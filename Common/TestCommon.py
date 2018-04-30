
import platform

WINDOWS = (platform.system() == 'Windows')

import Basic
import Program

Basic.onlyForInternalUse ()

print 'Test module Program.py'
Program.testAll ()

print '**** End of Common.py ****'