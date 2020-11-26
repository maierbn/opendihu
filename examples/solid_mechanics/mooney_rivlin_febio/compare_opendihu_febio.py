#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Script that parses the output from the febio and opendihu scenarios and computes a total error

import sys, os
import py_reader
import numpy as np

#print(os.getcwd())
# load data
data1 = py_reader.load_data(["build_release/out/febio_0000001.py"])
data2 = py_reader.load_data(["build_release/out/opendihu_0000001.py"])

# if files do not exist
if data1 == [] or data2 == []:
  quit()
  
component_name = "0"
total_error = 0

values1 = py_reader.get_values(data1[0], "geometry", component_name)
values2 = py_reader.get_values(data2[0], "geometry", component_name)


min_u1,max_u1 = py_reader.get_min_max(data1, "u", "2")
min_u2,max_u2 = py_reader.get_min_max(data2, "u", "z")

print("maximum z-displacement febio:    {}".format(max_u1))
print("maximum z-displacement opendihu: {}".format(max_u2))
  
# values2 contains entries for quadratic elements
# extract the corner values
n_elements = data2[0]['nElements']
nx = n_elements[0]
ny = n_elements[1]
nz = n_elements[2]
mx = nx*2 + 1
my = ny*2 + 1
mz = nz*2 + 1

values2_linear = []
for k in range(nz+1):
  for j in range(ny+1):
    for i in range(nx+1):
      values2_linear.append(values2[2*k * mx * my + 2*j * mx + 2*i])
  
#print("values1 (febio):    ",list(values1))
#print("values2 (opendihu): ",values2_linear)
  
error_rms = np.sqrt(np.mean((values1-values2_linear)**2))
  
print("rms: {}".format(error_rms))
