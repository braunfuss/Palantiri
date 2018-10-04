
import sys

import Basic
from   Program  import startTest                     # from Common

Basic.onlyForInternalUse()

sys.argv = startTest('cluster', workingDir='tools')

import cluster2                                      # from Cluster
cluster2.MainProc()

