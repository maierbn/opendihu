#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script reads a pair of *.exnode, *.exelem files and outputs a *.vtr file to be viewed with Paraview.
#
# usage: ./exfile_to_vtk.py <fibers.bin input file> <output_file (without .vtr)>

import sys, os
import numpy as np
import struct
import stl
from stl import mesh
import datetime
import pickle
import exnode_reader

from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle, VtkQuad
from pyevtk.hl import gridToVTK     # pip3 install pyevtk


input_filename_exnode = "Time_2_100.part0.exnode"
if len(sys.argv) >= 2:
  input_filename_exnode = sys.argv[1]
  
  input_filename_exelem,_ = os.path.splitext(input_filename_exnode)
  input_filename_exelem += ".exelem"
  
  if len(sys.argv) >= 3:
    input_filename_exelem = sys.argv[2]
else:
  print("usage: ./exfile_to_vtk.py <filename.exnode> [<filename.exelem>]")
  sys.exit(0)

output_filename,_ = os.path.splitext(input_filename_exnode)

print("filename exnode: \"{}\"".format(input_filename_exnode))
print("filename exelem: \"{}\"".format(input_filename_exelem))
print("output filename: \"{}.vtu\"".format(output_filename))

nodal_values = exnode_reader.parse_file(input_filename_exnode, [["coordinates",1],["coordinates",2],["coordinates",3]])

element_nodes = exnode_reader.parse_exelem_file(input_filename_exelem)

if True:
  print("element_nodes:")
  print(element_nodes)

  print("result:")
  print(nodal_values)

# create geometry for vtk unstructured grid

n_points = np.size(nodal_values,0)
n_elements = len(list(element_nodes))
print("n_points: {}, n_elements: {}".format(n_points, n_elements))

positions_x = np.array(nodal_values[:,0])
positions_y = np.array(nodal_values[:,1])
positions_z = np.array(nodal_values[:,2])

connectivity = np.zeros(n_elements*8)

for (element_no,element) in element_nodes.items():
  connectivity[(element_no-1)*8 + 0] = element[0]-1
  connectivity[(element_no-1)*8 + 1] = element[1]-1
  connectivity[(element_no-1)*8 + 2] = element[3]-1
  connectivity[(element_no-1)*8 + 3] = element[2]-1
  connectivity[(element_no-1)*8 + 4] = element[4]-1
  connectivity[(element_no-1)*8 + 5] = element[5]-1
  connectivity[(element_no-1)*8 + 6] = element[7]-1
  connectivity[(element_no-1)*8 + 7] = element[6]-1

offsets = np.array([(i+1)*8 for i in range(n_elements)])

cell_types = np.array([12 for _ in range(n_elements)])
element_nos = np.array(element_nodes.keys())

if False:
  print("connectivity: {}".format(connectivity))
  print("offsets: {}".format(offsets))
  print("cell_types: {}".format(cell_types))
  print("element_nos: {}".format(element_nos))

# write vtk file
unstructuredGridToVTK(output_filename, positions_x, positions_y, positions_z, connectivity = connectivity, offsets = offsets, \
                      cell_types = cell_types, cellData = {"element_no": element_nos}, pointData = None)
