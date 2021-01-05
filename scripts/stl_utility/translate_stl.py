#!/usr/bin/env python
# Translate mesh in stl file.

import numpy as np
import os, sys
import stl
from stl import mesh

# parse command line arguments
outfile = "out.stl"
infile = "in.stl"

if len(sys.argv) < 2:
  print("usage: {} <fileNameIn> <fileNameOut> <x> <y> <z>\n Translate vertices by [x,y,z].".format(sys.argv[0]))
  sys.exit(0)

if len(sys.argv) >= 2:
  if os.path.isfile(sys.argv[1]):
    infile = sys.argv[1]
  else:
    print("File \"{}\" does not exist.".format(sys.argv[1]))
    sys.exit(0)
  
if len(sys.argv) >= 3:
  outfile = sys.argv[2]
else:
  outfile = os.path.splitext(infile)[0]+"_out.stl"

# parse scaling factors
translation_x = 0.0
translation_y = 0.0
translation_z = 0.0

if len(sys.argv) == 6:
  translation_x = float(sys.argv[3])
  translation_y = float(sys.argv[4])
  translation_z = float(sys.argv[5])

translation = np.array([translation_x, translation_y, translation_z])
print("Input file:      \"{}\"".format(infile))
print("Output file:     \"{}\"".format(outfile))
print("Translation vector: [{}, {}, {}]".format(translation_x, translation_y, translation_z))
  
stl_mesh = mesh.Mesh.from_file(infile)

out_triangles = []

# loop over triangles in mesh
for p in stl_mesh.points:
  # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

  p1 = np.array(p[0:3])
  p2 = np.array(p[3:6])
  p3 = np.array(p[6:9])
  
  p1 += translation
  p2 += translation
  p3 += translation
  
  out_triangles += [[p1, p2, p3]]
  
n_triangles = len(out_triangles)

# Create the mesh
out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
for i, f in enumerate(out_triangles):
  out_mesh.vectors[i] = f
    
out_mesh.update_normals()
out_mesh.save(outfile)

print("File \"{}\" written.".format(outfile))
