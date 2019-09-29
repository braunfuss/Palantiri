import os
import sys


def start():

    if sys.argv[1] != 'new_version':
        return False

    at = os.path.join(os.getcwd(), 'Common', 'NewVersion.py')

    os.system(sys.executable + ' ' + at)
    return True
