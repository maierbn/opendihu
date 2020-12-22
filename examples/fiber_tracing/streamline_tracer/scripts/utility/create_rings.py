#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# This script extracts horizontal rings of edges from the stl mesh of biceps. The rings are not planar (but almost).
# Use create_rings.py instead, which samples the mesh at prescribed z values and produces planar rings.
# Output of this script is a pickle file `rings_created` that can be read in by ./create_mesh.py
#
# usage: ./create_rings.py <input file> [<n rings> [<min z> <max z>]]

import sys,os
import pickle
# add to python path in order to load stl_create_rings
if os.environ.get("OPENDIHU_HOME"):
  sys.path.append(os.path.join(os.environ.get("OPENDIHU_HOME"), "scripts/geometry_manipulation"))
import stl_create_rings

import datetime
now = datetime.datetime.now()
print(" ======= create_rings.py =======") 
print(now.strftime("%d/%m/%Y %H:%M:%S"))

# parse command line arguments
input_filename = "in.stl"

# parameters
# for initial mesh
bottom_clip = -600.0     # bottom clip plane of triangles that will not be considered
top_clip = -290.0        # top clip plane

# for repaired mesh ("biceps_full.stl")
bottom_clip = 37.0
top_clip = 300.0

n_loops = 20   # number of rings to extract

if len(sys.argv) < 2:
  print("usage: ./create_rings.py <input file> [<n rings> [<min z> <max z>]]")
  sys.exit(0)

if len(sys.argv) >= 2:
  if os.path.isfile(sys.argv[1]):
    input_filename = sys.argv[1]
  else:
    print("File \"{}\" does not exists".format(sys.argv[1]))
    sys.exit(0)
  
if len(sys.argv) >= 3:
  n_loops = int(sys.argv[2])
  print("n loops: {}".format(n_loops))
  
if len(sys.argv) >= 5:
  bottom_clip = float(sys.argv[3])
  top_clip = float(sys.argv[4])
  print("z range of rings: [{},{}]".format(bottom_clip, top_clip))
  
print("Input file: \"{}\"".format(input_filename))

write_output_mesh = True
loops = stl_create_rings.create_rings(input_filename, bottom_clip, top_clip, n_loops, write_output_mesh)

with open('rings_created', 'wb') as f:
  pickle.dump(loops, f)
  
