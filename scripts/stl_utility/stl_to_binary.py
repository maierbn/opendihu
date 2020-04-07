#!/usr/bin/env python

# @file stl_to_binary.py
# 
# Save the stl file in binary format
# 
# 
# @author Benjamin Maier
# 
# @date 01/2020
# 

import sys,os
import stl
from stl import mesh

# read in arguments
if len(sys.argv) < 2:
  print("usage: {} <input file> [<output file>]\n Show info.".format(sys.argv[0]))
  sys.exit(0)

# parse command line arguments
outfile = "out.stl"
infile = "in.stl"

if len(sys.argv) >= 2:
  if os.path.isfile(sys.argv[1]):
    infile = sys.argv[1]
  else:
    print("File \"{}\" does not exist.".format(sys.argv[1]))
    sys.exit(0)
  
if len(sys.argv) >= 3:
  outfile = sys.argv[2]
else:
  outfile = os.path.splitext(infile)[0]+"_binary.stl"

print("Input file: \"{}\"".format(infile))
print("Output file: \"{}\"".format(outfile))

stl_mesh = mesh.Mesh.from_file(infile)
stl_mesh.save(outfile)
