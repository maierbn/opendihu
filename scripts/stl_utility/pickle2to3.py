#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# To be executed with python3. Reads in a python2 pickle and outputs it as python 3.
#

import sys,os
import pickle

if len(sys.argv) < 3:
  print("usage: pickle2to3.py <input file> <output file>")
  sys.exit(0)

if len(sys.argv) >= 3:
  if os.path.isfile(sys.argv[1]):
    infile = sys.argv[1]
  else:
    print("File \"{}\" does not exist").format(sys.argv[1])
    sys.exit(0)
  
  outfile = sys.argv[2]
  
print("input file: \"{}\"".format(infile))
print("output file: \"{}\"".format(outfile))


# read in data
with open(infile, 'rb') as f:
  data = pickle.load(f, encoding='latin1')

# print data
print(data)

# output data
with open(infile, 'wb') as f:
  pickle.dump(data, f)
