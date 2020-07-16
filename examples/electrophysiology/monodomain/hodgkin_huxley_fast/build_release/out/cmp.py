#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import py_reader
import numpy as np

print(os.getcwd())
directory1 = "fast"
directory2 = "not_fast"

# read all files in directories
files1 = []
for filename in os.listdir(directory1):
  if filename.endswith(".py"):
    files1.append(os.path.join(directory1, filename))
files1 = sorted(files1)

files2 = []
for filename in os.listdir(directory2):
  if filename.endswith(".py"):
    files2.append(os.path.join(directory2, filename))
files2 = sorted(files2)

print("files: ",files1,files2)

# load data
data1 = py_reader.load_data(files1)
data2 = py_reader.load_data(files2)

n_values = min(len(data1), len(data2))

if len(data1) != len(data2):
  print("Warning: Directory {} contains {} files, directory {} contains {} files.".format(directory1, len(data1), directory2, len(data2)))

component_name = "0"
total_error = 0
for i in range(n_values):
  values1 = py_reader.get_values(data1[i], "solution", component_name)
  values2 = py_reader.get_values(data2[i], "solution", component_name)
  
  error = np.linalg.norm(values1-values2) / np.size(values1);
  total_error += error
  
  print("file no. {}, error: {}".format(i, error))
  
total_error /= n_values
print("avg error: {}".format(total_error))
