import os
import sys
from optparse import OptionParser

from palantiri.tools import config
from palantiri.common import Globals

parser = OptionParser(usage="%prog -f eventpath ")
parser.add_option("-f", "--evpath", type="string", dest="evpath", help="evpath")

(options, args) = parser.parse_args()


def init():

    C = config.Config(options.evpath)
    Config = C.parseConfig('config')

    tests = int(Config['runs'])
    import palantiri
    path = palantiri.__path__
    at = os.path.join(path[0], 'cluster/cluster2.py')
    cmd = sys.executable + ' ' + at + ' -f '+ options.evpath
    print('cmd = ', cmd)

    for i in range(tests):
        print('RUN: ', i)
        os.system(cmd)

    cmd = ('%s evaluatecluster.py -f %s') % (sys.executable,
                                             os.path.join(options.evpath,
                                                          'cluster'))

    at = os.path.join(path[0], 'cluster/evaluateCluster.py')
    cmd = sys.executable + ' ' + at + ' -f '+ os.path.join(options.evpath, "cluster")
    os.system(cmd)



def main():
    init()
