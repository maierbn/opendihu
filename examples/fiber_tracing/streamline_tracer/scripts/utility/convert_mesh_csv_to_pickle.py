#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# This script reads a csv file where in each line there are the node_positions of a fibre.
# It outputs the data as [[[n1x,n1y,n1z],[n2x,n2y,n2z],...], [...], ...].
#
# usage: ./convert_mesh_csv_to_pickle.py <input csv file> [<output file>]

import datetime
now = datetime.datetime.now()
print(" ======= convert_mesh_csv_to_pickle.py =======") 
print(now.strftime("%d/%m/%Y %H:%M:%S"))

import sys, os
import numpy as np
import pickle

infile = "streamlines.csv"
outfile = "out"
  
if len(sys.argv) < 2:
  print("usage: ./convert_mesh_csv_to_pickle.py <input csv file> [<output file>]")
  sys.exit(0)

if len(sys.argv) >= 2:
  if os.path.isfile(sys.argv[1]):
    infile = sys.argv[1]
  else:
    print("File \"{}\" does not exist.".format(sys.argv[1]))
    sys.exit(0)

if len(sys.argv) >= 3:
  outfile = sys.argv[2]
    
print("input file: {}".format(infile))
print("output file: {}".format(outfile))

# parse input file 
streamlines = []
with open(infile,"r") as f:
  for line in f:
    dataset = line.split(";")
    if dataset[-1] == "\n":
      dataset = dataset[:-1]
    dataset = map(float, dataset)
    
    streamline = []
    point = []
    for entry in dataset:
      point.append(entry)
      if len(point) == 3:
        streamline.append(point)
        point = []
        
    # compute length of streamline 
    length = 0
    p0 = None
    for point in streamline: 
      if p0 is not None:
        length += np.linalg.norm(np.array(p0)-np.array(point))
      p0 = point
    print("length: ", length)
    
    streamlines.append(streamline)
    
# dump data 
with open(outfile, 'wb') as f:
  pickle.dump(streamlines,f)
  
print("saved {} streamlines to {}".format(len(streamlines),outfile))
