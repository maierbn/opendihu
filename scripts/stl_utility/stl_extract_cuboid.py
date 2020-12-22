#!/usr/bin/env python
# Extract a cuboid portion from the mesh in the stl file.

import numpy as np
import os, sys
import stl
from stl import mesh

# parse command line arguments
outfile = "out.stl"
infile = "in.stl"

if len(sys.argv) < 2:
  print("usage: {} <fileNameIn> <fileNameOut> <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> [<mode>]".format(sys.argv[0]))
  print("  Extract mesh [xmin,xmax] x [ymin,ymax] x [zmin,zmax].")
  print("  mode: 0, 1 or 2, only triangles are extracted that have: 0=the center point, 1=any point, 2=all points inside")
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

xmin = -np.inf
xmax = np.inf
ymin = -np.inf
ymax = np.inf
zmin = -np.inf
zmax = np.inf

if len(sys.argv) >= 9:
  xmin = (float)(sys.argv[3])
  xmax = (float)(sys.argv[4])
  ymin = (float)(sys.argv[5])
  ymax = (float)(sys.argv[6])
  zmin = (float)(sys.argv[7])
  zmax = (float)(sys.argv[8])

mode = 0
if len(sys.argv) >= 10:
  mode = (int)(sys.argv[9])

print("Input file:      \"{}\"".format(infile))
print("Output file:     \"{}\"".format(outfile))
print("Extract mesh [{}, {}] x [{}, {}] x [{}, {}].".format(xmin,xmax,ymin,ymax,zmin,zmax))
print("mode: {} (only triangles are extracted that have: 0=the center point, 1=any point, 2=all points inside)".format(mode))

# load mesh  
stl_mesh = mesh.Mesh.from_file(infile)

out_triangles = []

# loop over triangles in mesh
for p in stl_mesh.points:
  # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

  p1 = np.array(p[0:3])
  p2 = np.array(p[3:6])
  p3 = np.array(p[6:9])
  
  bounding_box_xmin = min(p1[0],p2[0],p3[0])
  bounding_box_xmax = max(p1[0],p2[0],p3[0])
  bounding_box_ymin = min(p1[1],p2[1],p3[1])
  bounding_box_ymax = max(p1[1],p2[1],p3[1])
  bounding_box_zmin = min(p1[2],p2[2],p3[2])
  bounding_box_zmax = max(p1[2],p2[2],p3[2])
  center = 0.5*(p1+p2+p3)
  
  if mode == 0 and xmin <= center[0] <= xmax and ymin <= center[1] <= ymax and zmin <= center[2] <= zmax:
    out_triangles += [[p1, p2, p3]]
  elif mode == 1 \
    and ((xmin <= p1[0] <= xmax and ymin <= p1[1] <= ymax and zmin <= p1[2] <= zmax) or (xmin <= p2[0] <= xmax and ymin <= p2[1] <= ymax and zmin <= p2[2] <= zmax) or(xmin <= p3[0] <= xmax and ymin <= p3[1] <= ymax and zmin <= p3[2] <= zmax)):
    out_triangles += [[p1, p2, p3]]
  elif mode == 2 and xmin <= bounding_box_xmin and bounding_box_xmax <= xmax \
    and ymin <= bounding_box_ymin and bounding_box_ymax <= ymax \
    and zmin <= bounding_box_zmin and bounding_box_zmax <= zmax:
    out_triangles += [[p1, p2, p3]]
  
n_triangles = len(out_triangles)

# Create the mesh
out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
for i, f in enumerate(out_triangles):
  out_mesh.vectors[i] = f
    
out_mesh.update_normals()
out_mesh.save(outfile)

print("File \"{}\" written.".format(outfile))
