import os
import sys
from palantiri.common import Basic, Globals, Logfile


class MainObj(object):

    def __init__(self, externClass, version, runTimeLog=None, errorLog=None):

        self.version = version
        self.runTimeLog = runTimeLog
        self.errorLog = errorLog
        self.extern = externClass

    def run(self):

        Logfile.init(self.runTimeLog, self.errorLog)
        Logfile.setStartMsg(self.version)

        if not self.extern.init():
            Logfile.abort()

        try:
            ret = self.extern.process()

            if ret:
                msg = 'Palantiri finished'
            else:
                msg = 'Palantiri finished with error - maybe Sauron looked back?'

        except KeyboardInterrupt:
            msg = 'Gandalf made Pippin drop the Palantiri by Control C'
            ret = False

        self.extern.finish()
        Logfile.showLabel(msg)
        return ret
