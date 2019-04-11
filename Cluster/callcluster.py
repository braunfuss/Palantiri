import os
import sys
sys.path.append('../tools/')
sys.path.append('../Common/')

from optparse import OptionParser

import config
import Globals

parser = OptionParser(usage="%prog -f eventpath ")
parser.add_option("-f", "--evpath", type="string", dest="evpath", help="evpath")

(options, args) = parser.parse_args()


def init():

    C = config.Config(options.evpath)
    Config = C.parseConfig('config')

    tests = int(Config['runs'])

    at = os.path.join(os.getcwd(), 'cluster2.py')
    cmd = sys.executable + ' ' + at + ' -f '+ options.evpath
    print('cmd = ', cmd)

    for i in range(tests):
        print('RUN: ', i)
        os.system(cmd)

    cmd = ('%s evaluateCluster.py -f %s') % (sys.executable,
                                             os.path.join(options.evpath,
                                                          'cluster'))

    print(cmd)
    print(os.getcwd())
    os.system(cmd)


if __name__ == "__main__":

    init()
