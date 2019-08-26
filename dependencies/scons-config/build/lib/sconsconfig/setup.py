#!/usr/bin/env python
import subprocess, os, shutil

proj = os.path.basename(os.getcwd())
cmd = 'svn export config/proj %s'%proj
subp = subprocess.Popen(cmd, shell=True)
subp.wait()

for file in ['SConstruct', 'scons.py']:
    shutil.copy(os.path.join('config', file), '.')
