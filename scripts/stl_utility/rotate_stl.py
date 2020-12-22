#!/usr/bin/env python
# Rotate mesh in stl file.

import numpy as np
import os, sys
import stl
from stl import mesh

# parse command line arguments
outfile = "out.stl"
infile = "in.stl"

def print_usage():

  print("usage: \n\n"
    "  {p} <fileNameIn> <fileNameOut> <angle_x> <angle_y> <angle_z> [<x> <y> <z>]\n\n"
    "  {p} <fileNameIn> <fileNameOut> <axis_x> <axis_y> <axis_z> <angle> [<x> <y> <z>]\n\n"
    " The first version rotates by the given angles (in degrees) around the x,y and z-axis.\n"
    " The second version rotates by the given angle (in degree) around the vector [axis_x,axis_y,axis_z].\n\n"
    " The rotation center in both versions is the given center point [x,y,z] or the center of gravity if not specified.\n".format(p=sys.argv[0]))
  sys.exit(0)
  
# print usage
if len(sys.argv) < 2:
  print_usage()

# parse input filename
if len(sys.argv) >= 2:
  if os.path.isfile(sys.argv[1]):
    infile = sys.argv[1]
  else:
    print("File \"{}\" does not exist.".format(sys.argv[1]))
    sys.exit(0)
  
# parse output filename
if len(sys.argv) >= 3:
  outfile = sys.argv[2]
else:
  outfile = os.path.splitext(infile)[0]+"_out.stl"

# determine which of the two versions is used
version_no = 0
if len(sys.argv) == 6 or len(sys.argv) == 9:
  version_no = 1
elif len(sys.argv) == 7 or len(sys.argv) == 10:
  version_no = 2
else:
  print("Invalid number of arguments.\n")
  print_usage()

# find out center of gravity
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
  
center_of_gravity = points_sum / (3*len(stl_mesh.points))

print("Input file:         \"{}\"".format(infile))
print("Output file:        \"{}\"".format(outfile))
print("Center of gravity:   {}".format(center_of_gravity))

# version 1: rotate by the given angles (in degrees) around the x,y and z-axis
if version_no == 1:

  # parse angles
  x = float(sys.argv[3])
  y = float(sys.argv[4])
  z = float(sys.argv[5])
  angle_x = x/180.*np.pi
  angle_y = y/180.*np.pi
  angle_z = z/180.*np.pi

  rotation_x = np.array([[1, 0, 0], [0, np.cos(angle_x), -np.sin(angle_x)], [0, np.sin(angle_x), np.cos(angle_x)]])
  rotation_y = np.array([[np.cos(angle_y), 0, np.sin(angle_y)], [0, 1, 0], [-np.sin(angle_y), 0, np.cos(angle_y)]])
  rotation_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], [np.sin(angle_z), np.cos(angle_z), 0], [0, 0, 1]])
  
  rotation_matrix = rotation_z.dot(rotation_y.dot(rotation_x))

  # parse rotation point
  if len(sys.argv) == 9:
    rotation_point_x = float(sys.argv[6])
    rotation_point_y = float(sys.argv[7])
    rotation_point_z = float(sys.argv[8])
  else:      
    rotation_point_x = center_of_gravity[0]
    rotation_point_y = center_of_gravity[1]
    rotation_point_z = center_of_gravity[2]
      
  print("Rotations (degrees): {} {} {}".format(x, y, z))
  print("Rotation point:      [{}, {}, {}]".format(rotation_point_x, rotation_point_y, rotation_point_z))    
  
elif version_no == 2:
  # version 2: rotate by the given angle (in degree) around the vector [axis_x,axis_y,axis_z]
  
  # parse parameters
  axis_x = float(sys.argv[3])
  axis_y = float(sys.argv[4])
  axis_z = float(sys.argv[5])
  angle = float(sys.argv[6])
  angle = angle/180.*np.pi
  
  # normalize rotation axis
  axis = np.array([axis_x, axis_y, axis_z])
  axis = axis / np.linalg.norm(axis)
  axis_x = axis[0]
  axis_y = axis[1]
  axis_z = axis[2]
  
  # compute rotation matrix
  rotation_matrix = np.array([
    [np.cos(angle) + axis_x**2*(1 - np.cos(angle)), 
     axis_x*axis_y*(1 - np.cos(angle)) - axis_z*np.sin(angle), 
     axis_x*axis_z*(1 - np.cos(angle)) + axis_y*np.sin(angle)],
    [axis_y*axis_x*(1 - np.cos(angle)) + axis_z*np.sin(angle),
     np.cos(angle) + axis_y**2*(1 - np.cos(angle)),
     axis_y*axis_z*(1 - np.cos(angle)) - axis_x*np.sin(angle)],
    [axis_z*axis_x*(1 - np.cos(angle)) - axis_y*np.sin(angle),
     axis_z*axis_y*(1 - np.cos(angle)) + axis_x*np.sin(angle),
     np.cos(angle) + axis_z**2*(1 - np.cos(angle))]])
  
  # parse rotation point
  if len(sys.argv) == 10:
    rotation_point_x = float(sys.argv[7])
    rotation_point_y = float(sys.argv[8])
    rotation_point_z = float(sys.argv[9])
  else:      
    rotation_point_x = center_of_gravity[0]
    rotation_point_y = center_of_gravity[1]
    rotation_point_z = center_of_gravity[2]
  
  print("Rotation axis (normalized): [{}, {}, {}]".format(axis_x, axis_y, axis_z))
  print("Angle (degrees):     {}".format(angle*180./np.pi))
  print("Rotation point:      [{}, {}, {}]".format(rotation_point_x, rotation_point_y, rotation_point_z))  
  
rotation_point = np.array((rotation_point_x,rotation_point_y,rotation_point_z))
  
out_triangles = []
# rotate all points
# loop over triangles in mesh
for p in stl_mesh.points:
  # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

  p1 = np.array(p[0:3])
  p2 = np.array(p[3:6])
  p3 = np.array(p[6:9])
  
  # transform vertices
  vertex_list = []
  
  # apply rotation
  for vertex in [p1,p2,p3]:
    vertex = vertex - rotation_point
    vertex = rotation_matrix.dot(vertex)
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
