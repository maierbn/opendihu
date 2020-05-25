#!/usr/bin/env python
# Rotate mesh in stl file.

import numpy as np
import os, sys
import stl
from stl import mesh

# parse command line arguments
outfile = "out.stl"
infile = "in.stl"

if len(sys.argv) < 2:
  print("usage: {} <fileNameIn> <fileNameOut> <angle_x> <angle_y> <angle_z> [<x> <y> <z>]\n Rotate by angles (in degrees) by x,y,z-axis by [angle_x,angle_y,angle_z] around point [x,y,z] or center if not specified\n.".format(sys.argv[0]))
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

# find out center point
stl_mesh = mesh.Mesh.from_file(infile)
points_sum = np.array((0.0,0.0,0.0))
# loop over triangles in mesh
for p in stl_mesh.points:
  # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

  p1 = np.array(p[0:3])
  p2 = np.array(p[3:6])
  p3 = np.array(p[6:9])
  points_sum += p1
  points_sum += p2
  points_sum += p3
  
center = points_sum / (3*len(stl_mesh.points))

# parse scaling factors
rotation_point_x = center[0]
rotation_point_y = center[1]
rotation_point_z = center[2]

if len(sys.argv) >= 6:
  x = float(sys.argv[3])
  y = float(sys.argv[4])
  z = float(sys.argv[5])
  angle_x = x/180.*np.pi
  angle_y = y/180.*np.pi
  angle_z = z/180.*np.pi

if len(sys.argv) == 9:
  rotation_point_x = float(sys.argv[6])
  rotation_point_y = float(sys.argv[7])
  rotation_point_z = float(sys.argv[8])
  
rotation_point = np.array((rotation_point_x,rotation_point_y,rotation_point_z))
  
print("Input file:      \"{}\"".format(infile))
print("Output file:     \"{}\"".format(outfile))
print("Center:           {}".format(center))
print("Rotations (degrees): {} {} {}".format(x, y, z))
print("Rotation point:   {} {} {}".format(rotation_point_x, rotation_point_y, rotation_point_z))
  
out_triangles = []


# rotate all points
# loop over triangles in mesh
for p in stl_mesh.points:
  # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

  p1 = np.array(p[0:3])
  p2 = np.array(p[3:6])
  p3 = np.array(p[6:9])
  
  # transform vertices
  rotation_x = np.array([[1, 0, 0], [0, np.cos(angle_x), -np.sin(angle_x)], [0, np.sin(angle_x), np.cos(angle_x)]])
  rotation_y = np.array([[np.cos(angle_y), 0, np.sin(angle_y)], [0, 1, 0], [-np.sin(angle_y), 0, np.cos(angle_y)]])
  rotation_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], [np.sin(angle_z), np.cos(angle_z), 0], [0, 0, 1]])
  
  vertex_list = []
  
  for vertex in [p1,p2,p3]:
    vertex = vertex - rotation_point
    vertex = rotation_x.dot(vertex)
    vertex = rotation_y.dot(vertex)
    vertex = rotation_z.dot(vertex)
    vertex = vertex + rotation_point
    vertex_list.append(vertex)
  
  out_triangles += [vertex_list]
  
n_triangles = len(out_triangles)

# Create the mesh
out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
for i, f in enumerate(out_triangles):
  out_mesh.vectors[i] = f
    
out_mesh.update_normals()
out_mesh.save(outfile)

print("File \"{}\" written.".format(outfile))
