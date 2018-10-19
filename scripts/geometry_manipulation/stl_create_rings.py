#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# Sample horizontal rings/loops of the mesh at prescribed z values and produces planar rings.
# To use this function, call ./create_rings.py in system_tests fibers

import sys
import numpy as np
import csv
import collections
import copy
import scipy.spatial
import os
import pickle
import stl_create_mesh   # for standardize_loop

import stl
from stl import mesh
from svg.path import parse_path
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier

def get_intersecting_line_segment(triangle, z_value):
  """ return the line segment [pa, pb] of the triangle with z_value"""

  # barycentric coordinates
  # x(xi1,xi2) = (1-xi1-xi2)*x^{1} + xi1*x^{2} + xi2*x^{3},  xi1+xi2 <= 1, 0 <= xi1,xi2 <= 1
  # x_3(xi1,xi2) = z_value  =>  (1-xi1-xi2)*x^{1} + xi1*x^{2} + xi2*x^{3} = z_value
  #                         =>  xi1*(x^{2} - x^{1})  +  xi2*(x^{3} - x^{1})  =  z_value - x^{1}
  #                         =>  xi2 = ((z_value - x^{1}) - xi1*(x^{2} - x^{1})) / (x^{3} - x^{1})
  #                         =>  xi2 = (z_value - x^{1})/(x^{3} - x^{1}) - xi1 * (x^{2} - x^{1})/(x^{3} - x^{1})
  #                         =>  xi1 = (z_value - x^{1}) / (x^{2} - x^{1})  - xi2 * (x^{3} - x^{1}) / (x^{2} - x^{1}) 
  
  p1 = triangle[0]
  p2 = triangle[1]
  p3 = triangle[2]
  
  p1z = p1[2]
  p2z = p2[2]
  p3z = p3[2]

  # if p2z == p1z swap p1 and p3 
  if p2z == p1z:
    t = p1
    p1 = p3
    p3 = t
  
    p1z = p1[2]
    p3z = p3[2]
  
  debug = False
  has_to_intersect = False
  l = sorted([p1z, p2z, p3z, z_value])
  if l[0] != z_value and l[-1] != z_value:
    has_to_intersect = True
  #  debug = True
  
  if debug:
    print("points: {},{},{}, z_value: {}".format(p1z,p2z,p3z,z_value))
  
  # handle case where p1z = p2z
  if p2z == p1z:
    if p2z == z_value:
      return [p1, p2]
    if has_to_intersect:
      print("fail1, points: {},{},{}, z_value: {}".format(p1,p2,p3,z_value))
    return None
  
  # function that returns a point by its barycentric coordinates, xi1 and xi2
  def point(xi1,xi2):
    return (1-xi1-xi2) * p1 + xi1*p2 + xi2*p3
  
  # compute slope and xi1 offset of line segment in parameter space
  # xi1 = m*xi2 + c
  m = -(p3z - p1z) / (p2z - p1z)
  c = (z_value - p1z) / (p2z - p1z)
  
  if debug:
    print("xi1 = m*xi2 + c with m={}, c={}".format(m,c))
  
  # check which borders of the triangle in parameter space the line segment intersects
  intersects_xi2_equals_0 = (0 <= c <= 1)        # xi1 = c

  if abs(m) < 1e-12:
    intersects_xi1_equals_0 = False
  else:      
    intersects_xi1_equals_0 = (0 <= -c/m <= 1)     # m*xi2 + c = 0  =>  xi2 = -c/m
    
  if abs(1+m) < 1e-12:
    intersects_diagonal = False
  else:
    intersects_diagonal = (0 <= (c+m)/(1+m) <= 1)  # xi2 = 1-xi1  =>  xi1 = m*(1-xi1) + c  =>  c + m = (1+m)*xi1  =>  xi1 = (c+m)/(1+m)
    
  if debug:
    print("intersects_xi1_equals_0: {}, intersects_xi2_equals_0: {}, intersects_diagonal: {}, sum: {}".format(intersects_xi1_equals_0, intersects_xi2_equals_0, intersects_diagonal, int(intersects_xi2_equals_0) + int(intersects_xi1_equals_0) + int(intersects_diagonal)))
    
  # if triangle is not intersected by z plane
  if int(intersects_xi2_equals_0) + int(intersects_xi1_equals_0) + int(intersects_diagonal) != 2:
    
    if has_to_intersect:
      print("intersects_xi1_equals_0: {}, intersects_xi2_equals_0: {}, intersects_diagonal: {}, sum: {}".format(intersects_xi1_equals_0, intersects_xi2_equals_0, intersects_diagonal, int(intersects_xi2_equals_0) + int(intersects_xi1_equals_0) + int(intersects_diagonal)))
    
      print("fail2, points: {},{},{}, z_value: {}".format(p1z,p2z,p3z,z_value))
    return None
    
  # compute the intersection points
  edge = []
  if intersects_xi1_equals_0:
    xi1 = 0
    xi2 = -c/m
    edge.append(point(xi1,xi2))
    
  if intersects_xi2_equals_0:
    xi1 = c
    xi2 = 0
    edge.append(point(xi1,xi2))
    
  if intersects_diagonal:
    xi1 = (c+m) / (1+m)
    xi2 = 1 - xi1
    edge.append(point(xi1,xi2))
  
  return edge

def create_rings(input_filename, bottom_clip, top_clip, n_loops, write_output_mesh):
  """
  Create n_loops rings/loops (slices) on a closed surface, in equidistant z-values between bottom_clip and top_clip
  :param input_filename: file name of an stl file that contains the closed surface mesh of the muscle, aligned with the z-axis
  :param bottom_clip: the bottom z-value where to clip the muscle
  :param top_clip: the top z-value where to clip the muscle
  :param n_loops: number of loops/rings to extract between bottom_clip and top_clip
  :param write_output_mesh: if the result should be written to an output stl file for debugging
  """

  stl_mesh = mesh.Mesh.from_file(input_filename)

  out_triangles = []

  n_is_inside_1 = 0
  n_is_inside_2 = 0

  n_triangles = len(stl_mesh.points)

  z_samples = np.linspace(bottom_clip, top_clip, n_loops)

  # initialize list of empty lists for rings
  loops = []
  for i in range(n_loops):
    loops.append([])

  debug = False

  # loop over z samples
  for loop_no,z_value in enumerate(z_samples):
      
    print("loop no {}/{}".format(loop_no,len(z_samples)))

    # loop over all triangles in mesh
    for (no,p) in enumerate(stl_mesh.points):
      # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

      triangle = p

      p1 = np.array(triangle[0:3])
      p2 = np.array(triangle[3:6])
      p3 = np.array(triangle[6:9])
      
      # get intersection of current triangle with z plane
      edge = get_intersecting_line_segment([p1, p2, p3], z_value)
      
      # if there is an intersection
      if edge is not None:
      
        if debug:
          print("triangle ",p1,p2,p3,", found intersection edge: ", edge)
        
        # check if edge is already contained in loop
        edge_is_already_in_loop = False
        
        # for all points in the current loop
        for [p1, p2] in loops[loop_no]:
          if (np.allclose(p1, edge[0]) and np.allclose(p2, edge[1])) \
            or (np.allclose(p2, edge[0]) and np.allclose(p1, edge[1])):
            edge_is_already_in_loop = True
            break
        
        
        # append edge to loop
        if not edge_is_already_in_loop:
          if debug:
            print(" add edge ",edge," to loop no ", loop_no)
          #print(", prev: ", loops[loop_no], "->", 
          loops[loop_no].append(edge)
          #print(loops[loop_no]
        else:
          if debug:
            print("edge_is_already_in_loop")
        
  if debug:      
    print("-----------------")
    print("")
    print("loops: z values")
    for l in loops:
      print("")
      print("")
      print("")
      print(" -loop")
      print("")
      for e in l:
        print("  ",e)
      
  # sort edges of each loop
  # loop over z samples
  for loop_no in range(n_loops):
    
    loop = loops[loop_no]
    
    if len(loop) == 0:
      continue

    if debug:
      print("")
      print("loop ")
      print(loop[0:8])

    # find point with lowest x position as first_point
    first_edge = loop[0]
    first_point = first_edge[0]
    min_x = first_point[0]
    
    # iterate over points in current loop
    for edge in loop:
      for loop_point in edge:
        if loop_point[0] < min_x:
          min_x = loop_point[0]
          first_point = loop_point
      
    
    # fill new loop starting with first_point, ordered by consecutive edges
    new_loop = [first_point]
    
    if debug:
      print("first point: ", first_point)
    
    previous_end_point = first_point
    current_end_point = first_point
    
    #
    #new_loop = []
    #for edge in loop:
    #  for p in edge:
    #    new_loop.append(p)
    
    if True:
      while len(new_loop) < len(loop)+1:
          
        next_point_found = False
          
        # iterate over points in current loop
        for edge in loop:
          if np.allclose(current_end_point, edge[0]) and not np.allclose(previous_end_point, edge[1]):
            new_loop.append(edge[1])
            previous_end_point = current_end_point
            current_end_point = edge[1]
            next_point_found = True
            
            if debug:
              print("add point ",edge[1],", new length of loop: ",len(new_loop),", expected final length: ",len(loop)+1)
              
          elif np.allclose(current_end_point, edge[1]) and not np.allclose(previous_end_point, edge[0]):
            new_loop.append(edge[0])
            previous_end_point = current_end_point
            current_end_point = edge[0]
            next_point_found = True
            
            if debug:
              print("add point ",edge[0],", new length of loop: ",len(new_loop),", expected final length: ",len(loop)+1)
            
        # if the end point is again the start point, finish the loop. Note if by now still len(new_loop) < len(loop)+1 holds, there might be a different (distinct) loop for this z position which is discarded.
        if np.allclose(current_end_point, first_point):
          if debug:
            print("start point reached")
          break
          
        if not next_point_found:
          if debug: 
            print("no point found that continues loop")
          print("Error: loop for z={} could not be closed. Maybe there are triangles missing?".format(z_samples[loop_no]))
          break
            
    # store new loop 
    loops[loop_no] = list(new_loop)
        
  if debug:
    print("----------------")
    print("")
    print("loops: ",loops)
          
  print("The following rings have been extracted:")
  for (loop,z_value) in zip(loops,z_samples):
    print("at z = {}, n segments: {}".format(z_value,len(loop)))
          
  if write_output_mesh:
      
    # create output mesh
    out_triangles = []
    for loop in loops:
      #for edge in loop:
      #  p0 = edge[0]
      #  p1 = edge[1]
      #  p2 = 0.5*(p0+p1)
      #  out_triangles.append([p0, p1, p2])

      for i in range(len(loop)-1):
        p0 = loop[i]
        p1 = loop[i+1]
        p2 = 0.5*(p0+p1)
        out_triangles.append([p0, p1, p2])

      
    # create output mesh
    triangles = out_triangles
    n_triangles = len(out_triangles)
    print("n_triangles: ",n_triangles)

    # Create the mesh
    out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
    for i, f in enumerate(triangles):
      out_mesh.vectors[i] = f
      #for j in range(3):
        #print("set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
        
    #out_mesh.update_normals()

    outfile = "out/mesh_01.stl"
    out_mesh.save(outfile, mode=stl.Mode.ASCII)
    print("saved {} triangles to \"{}\" (loops)".format(n_triangles,outfile))

  # transform points from numpy array to list
  if False:
    list_loops = []
    for loop in loops:
      list_loop = []
      for point in loop:
        list_loop.append(point.tolist())
      list_loops.append(list_loop)
  
  return loops
