import os, subprocess
from command import *

def apply_patch(path, patch):
    cwd = os.getcwd()
    os.chdir(path)
    try:
        check_call('patch -N -p0 -i ' + patch)
    except:
        os.chdir(cwd)
        raise
    os.chdir(cwd)
