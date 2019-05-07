#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This scripts opens a binary file and interprets the contents as double values. It outputs the min/max and avg values.
#

import sys, os
import numpy as np
import struct

input_filename = "fibers.bin"

if len(sys.argv) == 1:
  print("usage: {} <filenames>".format(sys.argv[0]))
  sys.exit(0)
  
parameter_no = 1
while len(sys.argv) > parameter_no:
  input_filename = sys.argv[parameter_no]

  # read file
  values = []
  n_nans = 0
  n_infs = 0
  try:
    infile = open(input_filename, "rb")
    while True:
      double_raw = infile.read(8)
      if double_raw == b'' or len(double_raw) < 8:
        break
      value = struct.unpack('d', double_raw)[0]
      
      if np.isnan(value):
        n_nans += 1
      if np.isinf(value):
        n_infs += 1
        
      values.append(value)
      
    # output statistics
    print("file \"{}\", n_values: {}, n_nans: {}, n_infs: {}, min: {}, max: {}, mean: {}, first: {}, last: {}".\
      format(input_filename, len(values), n_nans, n_infs, min(values), max(values), np.mean(values), values[0], values[-1]))
  except:
    pass

  parameter_no += 1
