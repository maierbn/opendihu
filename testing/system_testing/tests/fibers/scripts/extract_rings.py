#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script extracts horizontal rings of edges from the stl mesh of biceps. The rings are not planar (but almost).
# Use create_rings.py instead, which samples the mesh at prescribed z values and produces planar rings.
#


import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import collections
import copy
import scipy.spatial
import os
import pymp
import pickle

import stl
from stl import mesh
from sets import Set
from svg.path import parse_path
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier

def get_edge_with_preferred_normal(p1, p2, p3, preferred_normal_direction):

  edge1 = -p1 + p2
  edge2 = -p2 + p3
  edge3 = -p3 + p1

  debug = False
  
  # compute normals on edges of triangle
  out1 = np.cross(edge3,edge1)
  normal1 = np.cross(out1,edge1)
  normal1 /= np.linalg.norm(normal1)
  
  out2 = np.cross(edge1,edge2)
  normal2 = np.cross(out2,edge2)
  normal2 /= np.linalg.norm(normal2)
    
  out3 = np.cross(edge2,edge3)
  normal3 = np.cross(out3,edge3)
  normal3 /= np.linalg.norm(normal3)
    
  # find out edge whose normal is the most in preferred_normal_direction
  dot1 = np.dot(normal1, preferred_normal_direction)
  dot2 = np.dot(normal2, preferred_normal_direction)
  dot3 = np.dot(normal3, preferred_normal_direction)
    
  dot_list = [abs(dot1),abs(dot2),abs(dot3)]
    
  if debug:
    print "preferred_normal_direction: ", preferred_normal_direction
    print "edge 1: {}, normal: {}, dot: {}".format(edge1, normal1, dot1)
    print "edge 2: {}, normal: {}, dot: {}".format(edge2, normal2, dot2)
    print "edge 3: {}, normal: {}, dot: {}".format(edge3, normal3, dot3)
    
  leading_edge = None
  if np.max(dot_list) == abs(dot1):
    if dot1 > 0:
      point0 = p1
      point1 = p2
    else:
      point0 = p2
      point1 = p1
  elif np.max(dot_list) == abs(dot2):
    if dot2 > 0:
      point0 = p2
      point1 = p3
    else:
      point0 = p3
      point1 = p2
  elif np.max(dot_list) == abs(dot3):
    if dot3 > 0:
      point0 = p3
      point1 = p1
    else:
      point0 = p1
      point1 = p3
  
  return point0,point1

# parse command line arguments
outfile = "out.stl"
infile = "in.stl"

if len(sys.argv) < 2:
  print "usage: extract_rings.py <input file> [<output file>]"
  sys.exit(0)

if len(sys.argv) >= 2:
  if os.path.isfile(sys.argv[1]):
    infile = sys.argv[1]
  else:
    print "File \"{}\" does not exists".format(sys.argv[1])
    sys.exit(0)
  
if len(sys.argv) >= 3:
  outfile = sys.argv[2]
else:
  outfile = os.path.splitext(infile)[0]+"_out.stl"

print "Input file: \"{}\"".format(infile)
print "Output file: \"{}\"".format(outfile)

stl_mesh = mesh.Mesh.from_file(infile)

out_triangles = []

n_is_inside_1 = 0
n_is_inside_2 = 0

n_triangles = len(stl_mesh.points)

# parameters
bottom_clip = -600.0     # bottom clip plane of triangles that will not be considered
top_clip = -290.0        # top clip plane

preferred_normal_direction = np.array([0.2,0.0,1.0])
preferred_normal_direction /= np.linalg.norm(preferred_normal_direction)

# extract "rings"/loops of edges from the triangle mesh
# the ring has a normal similar to preferred_normal_direction
loops = []

debug = False
  
# loop over all triangles in mesh as starting points for new rings
for (no,p) in enumerate(stl_mesh.points):
  # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

  print "no {}/{}".format(no,len(stl_mesh.points))

  triangle = p

  p1 = np.array(triangle[0:3])
  p2 = np.array(triangle[3:6])
  p3 = np.array(triangle[6:9])
  
  # exclude points that are to low (at the bottom end of the muscle)
  if p1[2] <= bottom_clip or p1[2] >= top_clip:
    continue
  
  if p1[2] < -575:
    preferred_normal_direction = np.array([0.2,0.0,1.0])
    preferred_normal_direction /= np.linalg.norm(preferred_normal_direction)
  else:
    preferred_normal_direction = np.array([0.0,0.0,1.0])
    preferred_normal_direction /= np.linalg.norm(preferred_normal_direction)


  
  # for debugging only use one no
  #if no != 1511:
  #if abs(p1[2] - -470.14199829) > 1:
  #  continue

  # get the egde [point1,point2] of the current triangle (p1,p2,p3) that has a normal most similar to preferred_normal_direction
  point1,point2 = get_edge_with_preferred_normal(p1, p2, p3, preferred_normal_direction)
  
  # check if point1 exists already in any other loop
  point_exists_in_loop = False
  for loop in loops:
    for point in loop:
      if np.allclose(point,point1):
        point_exists_in_loop = True
        break
    if point_exists_in_loop:
      break
  
  # if the considered point1 is already present in a loop, do not use it as starting point for a new loop (because it is not a new loop)
  if point_exists_in_loop:
    if debug:
      print "    point exists in loop"
    continue
  
  # now find other points of the "ring"/loop
  loop_points = []
  point_start = np.array(point1)    # copy starting point to later be able to find last point that closes the loop
  
  # store first two points of loop 
  loop_points.append(point_start)
  loop_points.append(point2)
  
  if debug:
    print "starting points of loop: ", point1, point2
  
  # loop until all points of the current loop are found
  end_of_loop = False
  while not end_of_loop:
    
    if debug:
      print "find triangle edge with first point ",point2
    triangle_found = False
    
    # loop again over all triangles, except the one from which the loop started
    for triangle0 in stl_mesh.points:
      # triangle contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

      # if this is the triangle where we started, do not consider it
      if np.allclose(triangle0, p):
        continue

      # extract nodes of triangle
      p1 = np.array(triangle0[0:3])
      p2 = np.array(triangle0[3:6])
      p3 = np.array(triangle0[6:9])
      
      # get the egde [point3,point4] of the current triangle (p1,p2,p3) that has a normal most similar to preferred_normal_direction
      # the searched edge is the one where point3 is the most recent point of the currently filled loop (the most recent point is stored in point2)
      point3,point4 = get_edge_with_preferred_normal(p1, p2, p3, preferred_normal_direction)
    
      if debug and False:
        print "  - check triangle with points ", point3,point4
    
      #print "point2: ",point2," check triangle with points3,4 ",point3,point4
    
      # check if point4 is already part of this loop. If it is, the current triangle is no longer considered ('continue' statement later)
      point4_exists_in_loop = False
      point3_exists_in_loop = False
      for point in loop_points:
        if np.allclose(point, point4):
          point4_exists_in_loop = True
        if np.allclose(point, point3):
          point3_exists_in_loop = True
          
      # check if the starting point, point3, of the found edge, [point3,point4], is the same as the most recent point of the loop, point2
      if np.allclose(point2, point3):
        if debug:
          print "yes found, point2 == point3, use next point4 ", point4, " point_start:",point_start
        
        if np.allclose(point4, point_start) and len(loop_points) == 2:
          if debug:
            print "avoid closing loop after 2 points"
          continue
        
        # if the other point is already in the current loop do not further consider this edge
        # if we did not reach the starting point but any other intermediate point of this loop, do not consider triangle edge, because it is a "short-cut" to the loop
        if point4_exists_in_loop and not np.allclose(point4,point_start):
          if debug:
            print "point exists already in loop and is not the start point"
          continue
        
        if debug:
          print "point ", point4
          
        # store new point to the loop
        loop_points.append(point4)
        
        # update most recent point
        point2 = point4
        
        # if we again reached the starting point, the loop is closed and we're done
        if np.allclose(point4, point_start):
          if debug:
            print "point4 == point_start"
          end_of_loop = True
          
        # we found the triangle that continues the loop, now successfully break the for loop over all triangles
        triangle_found = True
        break
        
      # check if the starting point, point3, of the found edge, [point3,point4], is the same as the most recent point of the loop, point2
      if np.allclose(point2, point4):
        if debug:
          print "yes found, point2 == point4, use next point3 ", point3, " point_start:",point_start

        if np.allclose(point3, point_start) and len(loop_points) == 2:
          if debug:
            print "avoid closing loop after 2 points"
          continue
          
        # if the other point is already in the current loop do not further consider this edge
        # if we did not reach the starting point but any other intermediate point of this loop, do not consider triangle edge, because it is a "short-cut" to the loop
        if point3_exists_in_loop and not np.allclose(point3,point_start):
          if debug:
            print "point exists already in loop and is not the start point"
          continue
          
        if debug:
          print "point ", point3    
              
        # store new point to the loop
        loop_points.append(point3)
        
        # update most recent point
        point2 = point3
        
        # if we again reached the starting point, the loop is closed and we're done
        if np.allclose(point3,point_start):
          if debug:
            print "point3 == point_start"
          end_of_loop = True
    
        # we found the triangle that continues the loop, now successfully break the for loop over all triangles
        triangle_found = True
        break
      
    # if the pass over all triangles did not result in a new point that continues the loop, end the search
    if not triangle_found:
      end_of_loop = True
      if debug:
        print "no further point found to close loop, point2: ", point2, ", point_start: ", point_start
    
  # store the currently collected loop to the list of loops
  if debug:
    print "    create loop with {} points: {}".format(len(loop_points), loop_points)
  loops.append(loop_points)
  
with open('rings_extracted', 'wb') as f:
  pickle.dump(loops, f)
  
# create output mesh
out_triangles = []
for loop in loops:
  for i in range(len(loop)-1):
    p0 = loop[i]
    p1 = loop[i+1]
    p2 = 0.5*(p0+p1)
    out_triangles.append([p0, p1, p2])
  
  
# create output mesh
triangles = out_triangles
n_triangles = len(out_triangles)
print "n_triangles: ",n_triangles

# Create the mesh
out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
for i, f in enumerate(triangles):
  out_mesh.vectors[i] = f
  #for j in range(3):
    #print "set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
    
#out_mesh.update_normals()

outfile = "mesh_01b.stl"
out_mesh.save(outfile, mode=stl.Mode.ASCII)
print "saved {} triangles to \"{}\" (loops)".format(n_triangles,outfile)
