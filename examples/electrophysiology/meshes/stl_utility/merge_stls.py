#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import struct
import stl
from stl import mesh
import shutil

# parse command line arguments
outfile = "out.stl"
infiles = []

if len(sys.argv) < 2:
  print "usage: show.py <input file1> <input file2> ... <output file>"
  sys.exit(0)

if len(sys.argv) >= 2:
  for arg in sys.argv[1:-1]:
    if os.path.isfile(arg):
      infiles += [arg]
    else:
      print "File \"{}\" does not exists".format(sys.argv[1])
  outfile = sys.argv[-1]
  if os.path.isfile(outfile):
    shutil.copyfile(outfile, outfile+'_backup')
  
print "Input files: \"{}\"".format(infiles)
print "Output file: \"{}\"".format(outfile)

result_mesh = None
for filename in infiles:
  
  new_mesh = mesh.Mesh.from_file(filename)
  if result_mesh is None:
    result_mesh = new_mesh
  else:
    result_mesh = mesh.Mesh(np.concatenate([result_mesh.data, new_mesh.data]))

result_mesh.save(outfile)
