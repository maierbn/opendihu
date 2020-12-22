#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# This script reads a csv file where in each line there are the node_positions of a fiber.
# It creates an STL file of the fibers that can be visualized, e.g. in paraview.

import datetime
now = datetime.datetime.now()
print(" ======= convert_mesh_csv_to_stl.py =======") 
print(now.strftime("%d/%m/%Y %H:%M:%S"))

import sys
import stl
from stl import mesh
import numpy as np


if len(sys.argv) >= 2:
  in_filename = sys.argv[1]
if len(sys.argv) >= 3:
  out_filename = sys.argv[2]
  
print("filename in: {}".format(in_filename))
print("filename out: {}".format(out_filename))

data = []
with open(in_filename,"r") as f:
  for line in f:
    dataset = line.split(";")
    if dataset[-1] == "\n":
      dataset = dataset[:-1]
    dataset = list(map(float, dataset))
    data.append(dataset)
#data = np.genfromtxt("streamlines.csv",delimiter=";")

#print data.shape

triangles = []

if type(data[0]) is not list:
  data = list([data])

for dataset in data:
  previous_point = None
  point = None
  #print "dataset --- "
  #print "dataset: ",dataset, len(dataset)
  for i in range(0,(len(dataset)//3)*3,3):
    previous_point = point
    point = np.array([dataset[i],dataset[i+1],dataset[i+2]])
    
    
    if previous_point is not None:
      #print("point ",[previous_point, point])
      triangles.append([previous_point, point, 0.5*(previous_point+point)])

# write debugging output stl meshes
def write_stl(triangles, outfile, description):
  # create output mesh
  n_triangles = len(triangles)

  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
    #for j in range(3):
      #print "set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
      
  out_mesh.update_normals()

  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\" ({})".format(n_triangles,outfile,description))

write_stl(triangles, out_filename, "streamlines")


