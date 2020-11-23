#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# This file contains functionality to creates a hexaeder mesh out of the rings given by create_rings.py. 
# The terms "ring" and "loop" are used synonymously in this file.
# Use create_mesh.py to actual perform the whole pipeline, or look at create_3d_mesh_simple or create_3d_mesh_from_border_points_faces.
# Output is a dict containing
#  {
#  "node_positions": nodes, 
#  "linear_elements": linear_elements, 
#  "quadratic_elements": quadratic_elements, 
#  "seed_points": seed_points,
#  "bottom_nodes": bottom_node_indices,
#  "top_nodes": top_node_indices,
#  "n_linear_elements_per_coordinate_direction": n_linear_elements_per_coordinate_direction,
#  "n_quadratic_elements_per_coordinate_direction": n_quadratic_elements_per_coordinate_direction,
#  }.
#
# usage: ./create_mesh.py [<triangulation_type> [<parametric_space_shape> [<n_points_x> [<n_grid_points_x>]]]]"

import sys, os
import numpy as np
import matplotlib

havedisplay = False
if not havedisplay:
  #print("use Agg backend")
  matplotlib.use('Agg')
else:
  print("use Tk backend")

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import collections, patches
import struct
import stl
from stl import mesh
import pickle
import scipy
import scipy.spatial
import scipy.linalg
import scipy.integrate
import scipy.optimize
import timeit

import stl_debug_output
import stl_create_rings
import spline_surface

def triangle_contains_point(triangle, point):
  """ check if a point lies inside a triangle, 2D """
  # barycentric coordinates
  # x(xi1,xi2) = (1-xi1-xi2)*x^{1} + xi1*x^{2} + xi2*x^{3},  xi1+xi2 <= 1, 0 <= xi1,xi2 <= 1
  # x_1(xi1,xi2) = point_1  =>  (1-xi1-xi2)*x^{1}_1 + xi1*x^{2}_1 + xi2*x^{3}_1 = point_1
  #                         =>  xi1*(x^{2}_1 - x^{1}_1)  +  xi2*(x^{3}_1 - x^{1}_1) =  point_1 - x^{1}_1
  # x_2(xi1,xi2) = point_2  =>  (1-xi1-xi2)*x^{1}_2 + xi1*x^{2}_2 + xi2*x^{3}_2 = point_2
  #                         =>  xi1*(x^{2}_2 - x^{1}_2)  +  xi2*(x^{3}_2 - x^{1}_2) =  point_2 - x^{1}_2

  # extract points as column vectors
  point0 = np.reshape(triangle[0],(-1,1))
  point1 = np.reshape(triangle[1],(-1,1))
  point2 = np.reshape(triangle[2],(-1,1))
  point = np.reshape(point,(-1,1))
  
  debug = False
  # compute xi1,xi2 of the point inside the triangle
  dudxi = np.concatenate([point1-point0, point2-point0],axis=1)
  try:
    dxidu = np.linalg.inv(dudxi)
  except:
    return (False, None)
  
  xi = dxidu.dot(point - point0)
  xi1 = xi[0]
  xi2 = xi[1]
  
  if debug:
    print("xi: {}, {}".format(xi1,xi2))
  
  # point is inside triangle iff conditions for xi1,xi2 are fulfilled
  tol = 1e-14
  condition = (-tol <= xi1 <= 1+tol and -tol <= xi2 <= 1+tol and xi1 + xi2 <= 1+tol)
  
  if debug:
    if condition:
      print("triangle: ({},{}) ({},{}) ({},{}), point: ({},{})".format(point0[0], point0[1], point1[0], point1[1], point2[0], point2[1], point[0], point[1]))
    
      print("dudxi: \n{}".format(dudxi))
      print("dxidu: \n{}".format(dxidu))
      print("p-p0: \n{}".format(point-point0))
      
      print("point found")
      print("")
  
  return (condition, xi)
  
def transform_to_world_space(x,y,triangles_parametric_space,triangle_list,parametric_space_shape):
  """
  map from the parametric domain to the world space
  :param x: first coordinate in the parametric domain, grid coordinate
  :param y: second coordinate in the parametric domain, grid coordinate
  :param triangles_parametric_space: the triangles in the parametric space that were mapped from the world space triangulation using the harmonic map
  :param triangle_list: list of triangles in the world space
  :param parametric_space_shape: which type of quadrangulation of the parameter space is given: 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid, 4 = like 3 but move grid points such that mean distance between points in world space gets optimal
  :return: point in world space that corresponds to (x,y) in parameter space
  """

  # transform to world space
  # find triangle in parametric space which contains grid point
  xi_point = None
  triangle_parameteric_space_no = None
  for (triangle_no,triangle_parametric_space) in enumerate(triangles_parametric_space):
    (contains_point, xi) = triangle_contains_point(triangle_parametric_space, np.array([x,y]))
    if contains_point:
      xi_point = xi
      triangle_parameteric_space_no = triangle_no
      break
      
  output_xi_msg = False
  if xi_point is None:
    for n_tries in range(5):
      if parametric_space_shape == 0 or parametric_space_shape == 3 or parametric_space_shape == 4:  # unit circle
        phi = np.arctan2(y,x)
        r = x / np.cos(phi)
        if (output_xi_msg): print(" (x,y) = ({},{}), (phi,r)=({}deg,{})  (check: (x,y)=({},{}))".format(x,y,phi*180./np.pi,r,r*np.cos(phi),r*np.sin(phi)))
        if abs(r) <= 1e-10:
          r_old = r
          r = 1e-4*np.sign(r)
          phi_old = phi
          phi += 1e-4
          x_old = x
          y_old = y
          x = r*np.cos(phi)
          y = r*np.sin(phi)
          if (output_xi_msg): print(" [1] adjust grid point (x,y) = ({},{})->({},{}), phi={}->{}, r={}->{}".format(x_old,y_old,x,y,phi_old,phi,r_old,r))
        
        elif abs(r) >= 1.0-1e-10:
          r_old = r
          r = 0.99*np.sign(r)
          phi_old = phi
          phi += 1e-4
          x_old = x
          y_old = y
          x = r*np.cos(phi)
          y = r*np.sin(phi)
          if (output_xi_msg): print(" [2] adjust grid point (x,y) = ({},{})->({},{}), phi={}->{}, r={}->{}".format(x_old,y_old,x,y,phi_old,phi,r_old,r))
        
        elif abs(r) <= 1e-4:
          r_old = r
          r = 0.01
          phi_old = phi
          phi = 0.0
          x_old = x
          y_old = y
          x = r*np.cos(phi)
          y = r*np.sin(phi)
          if (output_xi_msg): print(" [4] adjust grid point (x,y) = ({},{})->({},{}), phi={}->{}, r={}->{}".format(x_old,y_old,x,y,phi_old,phi,r_old,r))
        
        else:
          # move point a little
          r_old = r
          r = 0.99*r
          phi_old = phi
          phi += 1e-4
          x_old = x
          y_old = y
          x = r*np.cos(phi)
          y = r*np.sin(phi)
          if (output_xi_msg): print(" [3] adjust grid point (x,y) = ({},{})->({},{}), phi={}->{}, r={}->{}".format(x_old,y_old,x,y,phi_old,phi,r_old,r))
          # try again to find triangle in parametric space which contains grid point
          xi_point = None
          triangle_parameteric_space_no = None
          for (triangle_no,triangle_parametric_space) in enumerate(triangles_parametric_space):
            (contains_point, xi) = triangle_contains_point(triangle_parametric_space, np.array([x,y]))
            if contains_point:
              xi_point = xi
              triangle_parameteric_space_no = triangle_no
              break
      
          if xi_point is not None:
            if (output_xi_msg): print(" xi = [{},{}] found after {} tries".format(xi[0],xi[1],n_tries+1))
            break
    
  if xi_point is None:
    
    print("Error: could not find triangle in parameter space for grid point (x,y) = ({},{}), r={}".format(x,y,r))
    print("")
    return None
    
  if xi_point is not None:
      
    triangle_world_space = triangle_list[triangle_parameteric_space_no]
    p1 = triangle_world_space[0]
    p2 = triangle_world_space[1]
    p3 = triangle_world_space[2]
    
    xi = xi_point
    point_world_space = (1 - xi[0] - xi[1])*p1 + xi[0]*p2 + xi[1]*p3
    return point_world_space
  
def compute_mean_distances(grid_points_world_space,n_grid_points_x,n_grid_points_y):
  
  # evaluate quality of grid by computing mean distance between neighbouring points

  distances_current_loop = []
  relative_distances_current_loop = []
  
  # loop over grid points in world space
  for j in range(n_grid_points_y-1):
    for i in range(n_grid_points_x-1):
      
      # p6 p7 p8
      # p3 p4 p5
      # p0 p1 p2
      p4 = grid_points_world_space[j*n_grid_points_x+i]
        
      # p0
      if i > 0 and j > 0:
        p0 = grid_points_world_space[(j-1)*n_grid_points_x+i-1]
        distances_current_loop.append(np.linalg.norm(np.array(p4) - np.array(p0)))
            
      # p1
      if j > 0:
        p1 = grid_points_world_space[(j-1)*n_grid_points_x+i]
        distances_current_loop.append(np.linalg.norm(np.array(p4) - np.array(p1)))
        
      # p2
      if j > 0 and i < n_grid_points_x-1:
        p2 = grid_points_world_space[(j-1)*n_grid_points_x+i+1]
        distances_current_loop.append(np.linalg.norm(np.array(p4) - np.array(p2)))
      
      # p3
      if i > 0:
        p3 = grid_points_world_space[j*n_grid_points_x+i-1]
        distances_current_loop.append(np.linalg.norm(np.array(p4) - np.array(p3)))
      
      # p5
      if i < n_grid_points_x-1:
        p5 = grid_points_world_space[j*n_grid_points_x+i+1]
        distances_current_loop.append(np.linalg.norm(np.array(p4) - np.array(p5)))
      
      # p6
      if i > 0 and j < n_grid_points_y-1:
        p6 = grid_points_world_space[(j+1)*n_grid_points_x+i-1]
        distances_current_loop.append(np.linalg.norm(np.array(p4) - np.array(p6)))
      
      # p7
      if j < n_grid_points_y-1:
        p7 = grid_points_world_space[(j+1)*n_grid_points_x+i]
        distances_current_loop.append(np.linalg.norm(np.array(p4) - np.array(p7)))
      
      # p8
      if i < n_grid_points_x-1 and j < n_grid_points_y-1:
        p8 = grid_points_world_space[(j+1)*n_grid_points_x+i+1]
        distances_current_loop.append(np.linalg.norm(np.array(p4) - np.array(p8)))
  
      mean_distance_current_loop = np.mean(distances_current_loop)
      relative_distances_current_loop += [d / mean_distance_current_loop for d in distances_current_loop]
      #print("distances_current_loop: {}, mean: {}".format(distances_current_loop, mean_distance_current_loop)
      #print("relative_distances_current_loop: {}".format(relative_distances_current_loop)
      
  return distances_current_loop,relative_distances_current_loop

def standardize_loop(loop, lengths):
  """
  Rearrange the points in the loop such that the point with the lowest x value is the first and then the direction of the next points is counter-clockwise.
  :param loop: list of points, first and last point is the same
  :param lengths: the length of the loop gets computed and will be appended to the list lengths
  """
  
  debug = False
  
  # compute length of loop
  length = 0.0
  for i in range(len(loop)-1):
    length += np.linalg.norm(loop[i+1] - loop[i])
  lengths.append(length)

  # sort such that point with lowest x position is the first
  min_x = loop[0][0]
  start_point_index = 0
  for (no, [point_x, point_y, point_z]) in enumerate(loop):
    if point_x < min_x:
      min_x = point_x
      start_point_index = no
      
  if np.allclose(loop[-1], loop[0]):
    sorted_loop = loop[start_point_index:-1] + loop[0:start_point_index] + [loop[start_point_index]]
  else:
    sorted_loop = loop[start_point_index:] + loop[0:start_point_index]
  #sorted_loop = loop
  
  # find out orientation
  center_point = sum(sorted_loop)/len(sorted_loop)
  up = np.cross((-center_point + sorted_loop[0]),(-center_point + sorted_loop[1]))
  
  # if orientation is clockwise, reverse order of points in loop
  if np.sign(np.dot(up, np.array([0,0,1.0]))) < 0:
    sorted_loop = sorted_loop[::-1]
  
  # check if length is correct
  if debug:
    l = 0
    p0 = sorted_loop[0]
    for p in sorted_loop[1:]:
      l += np.linalg.norm(p-p0)
      p0 = p
    print("computed length: {}, stated length: {}".format(l, length))

  if debug:
    print("loop")
    print(loop)
    print("")
    
    print("sorted_loop")
    print(sorted_loop)
    print("")
    print("")
    
  return sorted_loop
    
def sample_border_points(loop, length, n_points, target_x, last_loop_start_point):
  """
  Given a loops of points (which are spaced like the initial STL triangles were), sample the loops at n_points equidistant points.
  The starting point is determined as follows: If last_loop_start_point is None, it determines the start point from target_x.
  The start point lies on the intersection of the plane x=target_x with the loop. if there are multiple intersections, the first that is in y- direction from the cog is taken.
  If last_loop_start_point is given, the start point is the point which has the lowest distance to last_loop_start_point., i.e. the connection line is perpendicular to the loop.
  This is helpful to ensure that two neighbouring loops on different z-levels yield similar border point loops.
  
  :param loop:   List of points that form the loop to be sampled.
  :param length:  The length of the loop, should be measured previously
  :param n_points: Number of points that will be created.
  :param target_x: a plane x=target_x that will determine the start point in the loop. The start point lies on that plane and is the one left (in x- direction) of the center.
  :param last_loop_start_point: If None, target_x will be used to determine the start_point, if given the start_point is chosen in a way such that its distance to last_loop_start_point is minimal.
  :return: start_point, border_points
  """
  
  debug = False
  border_points = []            # this list will hold the result
  h = float(length) / n_points  # parameter increment between border points
  
  # determine start point
  center_point = sum(loop)/len(loop)
  if target_x is None:
    target_x = 0.5*(loop[0][0]+loop[-1][0])   # plane x=target_x on which the start point will lie, it has to be ensured, that all rings touch this plane
  
  if debug:
    print("sample_border_points, loop with {} points, length: {}, sample {} points, target_x: {}, h: {}, last_loop_start_point:{}\n".format(len(loop), length, n_points, target_x, h, last_loop_start_point))
    print("loop: {}".format(loop))
    
    # check if given length is correct
    l = 0
    p0 = loop[0]
    for p in loop[1:]:
      l += np.linalg.norm(p-p0)
      p0 = p
    print("computed length: {}, stated length: {}".format(l, length))
  
  start_point = None
  start_point_point_index = None
  
  # determine if last_loop_start_point is used
  use_last_loop_start_point = False
  if last_loop_start_point is not None:
    use_last_loop_start_point = True
    
  # determine start point from last_loop_start_point
  if use_last_loop_start_point:
    last_loop_start_point = np.array(last_loop_start_point)
    
    # add first point as candidate in case none is found
    start_point_candidates = [[np.linalg.norm(np.array(loop[0])-last_loop_start_point), loop[0], 0]]
    
    # loop over loop points
    for point_index in range(1,len(loop)):
      p0 = np.array(loop[point_index-1])
      p1 = np.array(loop[point_index])
      u = -p0 + p1   # edge
      
      # determine loop_start_point as the point in the loop that is closest to the given last_loop_start_point
      if u.dot(u) < 1e-10:
        t_start = -1
      else:
        t_start = (last_loop_start_point - p0).dot(u) / u.dot(u)
      
      if debug:
        print("t: {}".format(t_start))
      
      if t_start >= 0 and t_start <= 1:
        plumb_foot_point_start = p0 + t_start * u
        
        distance = np.linalg.norm(plumb_foot_point_start-last_loop_start_point)
        start_point_candidates.append([distance, plumb_foot_point_start, point_index-1])
        
        if debug:
          print("at point_index {}, t: {}, plumb_foot_point_start: {}".format(point_index, t_start, plumb_foot_point_start))
      elif t_start < 0 and t_start != -1:
        point = p0
        
        distance = np.linalg.norm(point-last_loop_start_point)
        start_point_candidates.append([distance, point, point_index-1])
      elif t_start > 1:
        point = p1
        
        distance = np.linalg.norm(point-last_loop_start_point)
        start_point_candidates.append([distance, point, point_index-1])
        
    
    # determine candidate with lowest distance
    start_point_candidates.sort(key = lambda item: item[0])
    
    if debug:
      print("last_loop_start_point: {}, start_point_candidates (distance,point,index): {}".format(last_loop_start_point,start_point_candidates))
    
    # now we use the candidate with the lowest distance
    start_point_best_candidate = start_point_candidates[0]
    
    start_point = start_point_best_candidate[1]
    start_point_point_index = (start_point_best_candidate[2]+1)%len(loop)
    point_index = start_point_best_candidate[2]
    
  else:
    point_index = 0
    
  previous_loop_point = loop[point_index]
  point_index = (point_index+1)%len(loop)
  start_point_reached = False
  loop_iteration_no = 0
  end_iteration = False
  
  # iterate over the loop of points until the start point is reached, then iterate one more time, collecting the sampled border points in equidistant intervals
  while not end_iteration:
    loop_point = loop[point_index]
    
    if debug:
      print(" - loop iteration point_index: {}, start_point_point_index: {} (start_point: {})".format(point_index, start_point_point_index, start_point))
    
    # (if start_point was not determined) if the current edge contains a point that is horizontal in y-direction to the center point and left of it
    # or (if start_point was determined) if the start point is the first one
    if (not use_last_loop_start_point \
      and previous_loop_point[1] < center_point[1] and loop_point[1] < center_point[1] \
      and (previous_loop_point[0] <= target_x <= loop_point[0] or \
        loop_point[0] <= target_x <= previous_loop_point[0])) \
        \
      or (use_last_loop_start_point and (not start_point_reached or point_index == start_point_point_index)):
    #if previous_loop_point[1] < center_point[1] and loop_point[1] < center_point[1] \
    #  and (previous_loop_point[0] <= target_x <= loop_point[0] or \
    #    loop_point[0] <= target_x <= previous_loop_point[0]):
        
      # here, we are at the start point
        
      # if we meet this point again, the loop is closed, 
      # then proceed with current edge, p0-p1 where is ps is on it (p0-ps-p1), because the part p0-ps was not yet visited
      if start_point_reached:
        end_iteration = True
      else:
        # now the search for the start point has ended and the normal iterations where border points will be sampled continues
        start_point_reached = True
        
        # if the start point is not yet set
        if start_point is None:
          edge = -previous_loop_point + loop_point 
          edge_length = np.linalg.norm(edge)
          alpha = (target_x - previous_loop_point[0]) / edge_length
          start_point = previous_loop_point + alpha*edge
        
        # add starting point of loop as first border points
        border_points.append(start_point)
        if debug:
          print("   add start point: ",border_points)
        
        previous_loop_point = start_point
        
        t_previous_loop_point = 0
        t_next_border_point = h
          
    # if we are at the second part of the loop iteration where the actual border points are collected
    if start_point_reached:
      # compute current edge
      edge = -previous_loop_point + loop_point 
      edge_length = np.linalg.norm(edge) 
      
      n_on_edge = 0
      
      # collect all border points on current edge
      while t_next_border_point <= t_previous_loop_point + edge_length:
        
        border_point = previous_loop_point + edge * (t_next_border_point - t_previous_loop_point) / edge_length
        border_points.append(border_point)
        
        if debug:
          print("   > border point {}, t previous: {}, edge length: {}, next: {},  t_previous_loop_point + edge_length: {}".format(border_point, t_previous_loop_point, edge_length, t_next_border_point,  t_previous_loop_point + edge_length))
        
        t_next_border_point += h
        if debug:
          n_on_edge += 1
          print("   new t_next_border_point: ",t_next_border_point)
      
      if debug:
        print("   edge length {}, edge ({}, {}), n_on_edge: {}, t_previous_loop_point: {} t_next_border_point: {}".format(edge_length,previous_loop_point,loop_point,n_on_edge,t_previous_loop_point,t_next_border_point))
      
      # move on to next edge
      t_previous_loop_point += edge_length
      
    previous_loop_point = loop_point
    point_index = (point_index+1)%len(loop)
    
    loop_iteration_no += 1
    
    # if the start point was not found
    if loop_iteration_no == 2*len(loop):
      print("Warning: target plane y={} does not intersect loop, center_point: {}!".format(target_x, center_point))
      target_x = None
      loop_iteration_no = 0
      last_loop_start_point = loop[0]
      
      # restart with last_loop_start_point set to the first point in the loop
      return sample_border_points(loop, length, n_points, target_x, last_loop_start_point)
  
  # end of while loop
  
  if debug:
    print("visited ring length: {}, actual length: {}".format(t_previous_loop_point, length))
  
  if debug:
    print("distances between border points, note this are not necessary the same as h:")
    p0 = border_points[0]
    for p in border_points[1:]:
      print("  {}".format(np.linalg.norm(p0-p)))
      p0 = p
  
  # If there were too many points collected, crop to needed points.
  # This happens because on the edge of the starting point, p0-ps-p1 (with edge p0-p1 and ps=starting point), the part ps-p1 is visited twice.
  if len(border_points) > n_points:
    if debug:
      print("too many points: {}, n_points: {}, {}\n\n".format(len(border_points), n_points, border_points))
    border_points = border_points[:n_points]
  elif len(border_points) < n_points:
    print("Error, not enough points: {}, n_points: {}, {}\n\n".format(len(border_points), n_points, border_points))
  
  if debug:
    print("result: ")
    print("start point: ",start_point)
    print("border points: ",len(border_points))
    print(border_points)
    print("h: ",h)
    print("n border points: {}".format(len(border_points)))
  
  return start_point,border_points
    
def rings_to_border_points(loops, n_points):
  """
  Standardize every ring to be in counter-clockwise direction and starting with the point with lowest x coordinate,
  then sample border points
  :param loops: list of rings
  :param n_points: number of border points that should be sampled per ring
  :returns: border_point_loops,lengths  where lengths is a list of lengths of each loop and border_point_loops is the list of loops with border points
  """
  # standardize all rings
  sorted_loops = []
  lengths = []
  for loop in loops:

    if not loop:
      continue

    sorted_loop = standardize_loop(loop, lengths)
    sorted_loops.append(sorted_loop)
  
  # get center of gravity of the cogs of each loop
  #mean_point = sum(sum(loop)/len(loop) for loop in loops)/len(loops)
  #target_x = mean_point[1]
  
  # sample border points for each ring
  border_point_loops = []
  last_loop_start_point = None
  for loop_no,(loop,length) in enumerate(zip(sorted_loops,lengths)):
    
    if not loop:
      continue
    
    # determine target_x, the start point of the loop will lie on the plane x=target_x
    center_points = []
    smoothing_width = min(loop_no, len(sorted_loops)-loop_no-1, 5)
    for i in range(-smoothing_width,smoothing_width+1):
      loop_index = loop_no + i
      if 0 <= loop_index < len(sorted_loops):
        loop0 = sorted_loops[loop_index]
        center_point = sum(loop0)/len(loop0)
        center_points.append(center_point)
    
    center_point = sum(center_points)/len(center_points)
    target_x = center_point[0]
    
    # sample loop with n_points equidistant points
    start_point,border_points = sample_border_points(loop, length, n_points, target_x, last_loop_start_point)
    last_loop_start_point = start_point
    border_point_loops.append(border_points)
    
  return border_point_loops,lengths

def border_point_loops_to_list(border_point_loops):
  """
  transform the points from numpy array to list, such that they can be extracted from the opendihu C++ code
  """
  result = []
  for border_point_loop in border_point_loops:
    result_border_points = []
    for border_point in border_point_loop:
      result_border_points.append(border_point.tolist())
    result.append(result_border_points)
  return result

def create_planar_mesh(border_points, loop_no, n_points, \
    n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, improve_mesh, show_plot, debugging_stl_output, stl_triangle_lists):
  """
  Create a 2d mesh using the harmonic map algorithm. Input is a list of border points which will be the boundary nodes of the mesh.
  :param border_points: list of points on the border that will be the nodes, (each point is again a list with 3 entries)
  :param loop_no: the number of the current loop, this is only needed for debugging output filename
  :param n_points: number of points on the border of the resulting mesh
  :param triangulation_type: 0 = scipy, 1 = triangle, 2 = center pie (2 is best), 3 = minimized distance
  :param parametric_space_shape: 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
  :param max_area_factor: only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
  :param improve_mesh: if the generated mesh should be smoothed afterwards
  :param show_plot: if the matplotlib plots of geometry and harmonic map should be generated
  :param debugging_stl_output: if list should be filled with STL triangles that can be output to a STL mesh for debugging
  :param stl_triangle_lists: the debugging lists: [out_triangulation_world_space, markers_border_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,markers_grid_points_parametric_space, markers_grid_points_world_space]
  """

  # note: the first and last point in the loop is the same, because it represents the edges between the points

  global distances_between_world_mesh_nodes_std   # list of distances between neighboring grid points
  global relative_distances_between_world_mesh_nodes_std   # list of relative distances between neighbouring grid points, relative per slice
  
  t_start = timeit.default_timer()

  debug = False   # enable debugging output
  
  if debugging_stl_output:
    [out_triangulation_world_space, markers_border_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,\
      markers_grid_points_parametric_space, markers_grid_points_world_space] = stl_triangle_lists
    
  if len(border_points) != n_points:
    print("\nError in create_planar_mesh, n_points: {}, but {} border points given.\n border points: {}\n".format(n_points, len(border_points), border_points))
    
  # triangulate surface in world space
  points = np.reshape(border_points,(n_points,3))
  
  if debug:  
    print("")
    print("points:")
    print(points)

  if debugging_stl_output:
    # for debugging create markers at border points
    factor = 1.0
    for point in points:
      factor *= 1.02
      if factor >= 3:
        factor = 3
      size = 0.1*factor
      diag0 = np.array([-size,-size,-size])
      diag1 = np.array([size,-size,-size])
      diag2 = np.array([-size,size,-size])
      diag3 = np.array([size,size,-size])
      diag4 = np.array([-size,-size,size])
      diag5 = np.array([size,-size,size])
      diag6 = np.array([-size,size,size])
      diag7 = np.array([size,size,size])
      markers_border_points_world_space += [
        [point+diag0,point+diag3,point+diag1],[point+diag0,point+diag2,point+diag3],  # bottom
        [point+diag4,point+diag5,point+diag7],[point+diag4,point+diag7,point+diag6],  # top
        [point+diag0,point+diag1,point+diag5],[point+diag0,point+diag5,point+diag4],  # front
        [point+diag2,point+diag7,point+diag3],[point+diag2,point+diag6,point+diag7],  # back
        [point+diag2,point+diag0,point+diag4],[point+diag2,point+diag4,point+diag6],  # left
        [point+diag1,point+diag3,point+diag7],[point+diag1,point+diag7,point+diag5]  # right
      ]

  # extract z value 
  center_point = np.sum(points,axis=0)/len(points)
  z_value = center_point[2]
  (max_x,max_y,max_z) = np.max(points,axis=0)
  (min_x,min_y,min_z) = np.min(points,axis=0)
  extent_x = max_x - min_x
  extent_y = max_y - min_y
  
  # store points, because later they will be overwritten by adding new points from the triangulation
  original_points = np.array(points)
  n_original_points = n_points
  n_points_per_face = (int)(n_original_points/4)
  
  # project points on xy=z_value plane
  projected_points = []
  for point in points:
    projected_points.append(np.array([point[0], point[1]]))
  
  projected_points = np.reshape(projected_points, (-1,2))

  # define helper variables with default values that are only later used when triangulation_type == 2
  def get_modified_phi(phi_in): 
    phi_out = phi_in
    return phi_out
  modify_phi = False
  determine_additional_points_on_ring = False
  n_additional_points_on_ring = 0
  

  if debug:
    print("")
    print("projected_points:")
    print(projected_points)
  
  if triangulation_type == 0:
    # delaunay triangulation of scipy, does not add new points but is not constrained (works for convex domains)
  
    triangulation = scipy.spatial.Delaunay(projected_points)
    point_indices_list = triangulation.simplices
    triangle_list = points[point_indices_list]
    
  elif triangulation_type == 1:
    # delaunay triangulation of triangle package, adds new points, is constrained (works for concave domains)
    
    import triangle   # sudo pip triangle

    # create delaunay triangulation of points
    segments = np.reshape([[i,i+1] for i in range(n_points)], (n_points,2))
    segments[n_points-1] = np.array([n_points-1,0])
    
    data = {"vertices": projected_points, "segments": segments}
  
    max_area = extent_x * extent_y / max_area_factor
    if debug:
      print("maximum area: ",max_area)
  
    #triangulation = triangle.triangulate(data, 'pq')
    
    try:
      triangulation = triangle.triangulate(data, 'pqa'+str(max_area))
    except:
      print("Triangulation failed, n_points: {}, max_area: {}, loop_no: {}, projected_points: {}".format(n_points,max_area,loop_no,projected_points))
      
      t_stop = timeit.default_timer()
      duration = t_stop - t_start
      return [], duration
      
    triangulated_projected_points = np.array(triangulation['vertices'])
    
    # transform projected points back to 3D points
    points = []
    for projected_point in triangulated_projected_points:
      points.append(np.array([projected_point[0], projected_point[1], z_value]))
    
    # update n_points
    n_points = len(points)
    points = np.reshape(points, (-1,3))
    
    # determine the phi angle in the circle of the current point
    original_point_phi_value = []
    for original_point_no in range(n_original_points):
      border_point = points[original_point_no]
    
      phi = float(original_point_no) / n_original_points * 2 * np.pi
      original_point_phi_value.append(phi)
    
    # add additional points on ring
    # settings
    determine_additional_points_on_ring = True
    rescale_phi = True
    if parametric_space_shape == 3 or parametric_space_shape == 4:
      modify_phi = True
    
    # normal implementation without searching for additional border points on ring that the triangulation created
    if not determine_additional_points_on_ring:
      point_indices_list = triangulation["triangles"]
      triangle_list = points[point_indices_list]
      n_additional_points_on_ring = 0
      
    else:
      
      additional_points_on_ring = []
      new_points = list(points[0:n_original_points])
      interior_points = []
      
      # loop new points that were added by triangulation and are not the initial boundary points 
      for point_index in range(n_original_points,n_points):
        point = points[point_index]
        
        current_point_is_on_boundary = False
        # check if this point lies between two border points
        # loop over initial border points
        for border_point_index in range(n_original_points):
          border_point1 = points[border_point_index]
          border_point2 = points[(border_point_index+1)%n_original_points]
          
          v1 = -point + border_point1
          v2 = -point + border_point2
          v3 = -border_point1 + border_point2
          
          v1 = v1 / np.linalg.norm(v1)
          v2 = v2 / np.linalg.norm(v2)
          
          # if the point lies on the line between the two border points
          if abs(np.linalg.norm(np.cross(v1,v2))) < 1e-3:
            if abs(v3[0]) < abs(v3[1]):            
              alpha = (point[1] - border_point1[1]) / v3[1]
            else:
              alpha = (point[0] - border_point1[0]) / v3[0]
            
            if alpha > 1.0 or alpha < 0.0:
              #print("alpha: {} continue".format(alpha))
              continue
            
            phi = original_point_phi_value[border_point_index] + alpha * \
               (original_point_phi_value[(border_point_index+1)%n_original_points] - original_point_phi_value[border_point_index])
            original_point_phi_value.append(phi)
        
            #print("{} alpha: {}, phi: {} {} -> {}".format(point_index, alpha, original_point_phi_value[border_point_index], original_point_phi_value[(border_point_index+1)%n_original_points], phi))
        
            new_points.append(point)
            additional_points_on_ring.append(point_index)
            current_point_is_on_boundary = True
            break
        
        if not current_point_is_on_boundary:
          interior_points.append(point)
      
      # store points such that points = [<previous original points>, <newly determined points on the border>, <rest of points>]
      
      original_points = np.array(new_points)
      new_points += interior_points
      
      n_additional_points_on_ring = len(additional_points_on_ring)
      
      #print("n_additional_points_on_ring: {}".format(n_additional_points_on_ring))
      
      # adjust triangulation indices
      point_indices_list = triangulation["triangles"]
      
      for i in range(len(point_indices_list)):
        for point_no in range(len(point_indices_list[i])):
          point_index = point_indices_list[i][point_no]
          
          if point_index >= n_original_points:
            
            # count number of border points before old point_index
            n_additional_points_on_ring_before_point_index = 0
            for j in range(n_additional_points_on_ring):
              if additional_points_on_ring[j] < point_index:
                n_additional_points_on_ring_before_point_index += 1
              elif additional_points_on_ring[j] == point_index:
                point_indices_list[i][point_no] = n_original_points + n_additional_points_on_ring_before_point_index
                
                break
              else:
                point_indices_list[i][point_no] = point_index + n_additional_points_on_ring - n_additional_points_on_ring_before_point_index
                break
      
      # points has the following structure: [<list of original border points>, <list of new border points>, <list of interior points>]
      # original_points has the following structure: [<list of original border points>, <list of new border points>]
      points = np.array(new_points)
      triangle_list = points[point_indices_list]
      n_regular_grid_border_points = n_original_points
      n_original_points += n_additional_points_on_ring
      
      previous_original_point_phi_value = list(original_point_phi_value)
      
      # set phi values
      if rescale_phi:
        indices = np.argsort(original_point_phi_value)
        #print("original_point_phi_value: {}".format(original_point_phi_value))
        equidistant_values = np.linspace(0, 2*np.pi, n_original_points+1)[0:-1]
        #print("equidistant_values: {}".format(equidistant_values))
        #print("indices: {}".format(indices))
        for i,index in enumerate(indices):
          original_point_phi_value[index] = equidistant_values[i]
        #print("original_point_phi_value: {}".format(original_point_phi_value))
    
      #print("previous border points: {}, n_original_points: {}, n_additional_points_on_ring: {}, interior: {}, n_points: {}".\
      #  format(n_original_points-n_additional_points_on_ring, n_original_points, n_additional_points_on_ring, n_points-n_original_points, n_points))
      #print("additional_points_on_ring: {}".format(additional_points_on_ring))
      
      # setup map between parameter space regular grid in the circle and the transformed parameter space grid with the additional border points
      # this is done by defining a map for phi
      # map from phi to phi
      
      def get_modified_phi(phi_in):
        
        # normalize phi to [0,2*np.pi)
        if phi_in < 0:
          phi_in += 2*np.pi
        
        # determine position of phi between regular grid border points
        phi_increment = (2*np.pi) / n_regular_grid_border_points
        previous_border_point_index = (int)(phi_in / phi_increment)
        
        # determine factor between previous and next border point
        alpha = (phi_in - previous_border_point_index*phi_increment) / phi_increment
        
        # determine positions of phi in the new border points
        next_phi_value = 2*np.pi
        if previous_border_point_index+1 < len(original_point_phi_value):
          next_phi_value = original_point_phi_value[previous_border_point_index+1]
          
        previous_phi_value = original_point_phi_value[previous_border_point_index]
        
        # compute phi value with alpha between new border points
        phi_out = previous_phi_value + alpha * (next_phi_value - previous_phi_value)
        
        #print("phi_in: {}, phi_increment: {}, previous_border_point_index:{} [{},{}], alpha:{} new:[{},{}], phi_out: {}".format(phi_in, phi_increment, previous_border_point_index, previous_border_point_index*phi_increment, (previous_border_point_index+1)*phi_increment, alpha,\
        #  previous_phi_value, next_phi_value, phi_out))
        
        return phi_out
    
  elif triangulation_type == 2 or triangulation_type == 3:
    # 2: simple custom triangulation with triangles around one center point in CoG
    # 3: custom triangulation with triangles around point for which distance is minimized

    # compute the center point by minimizing the distances to the border points
    if triangulation_type == 3:
      
      # objective function
      def squared_distance_to_all_points(center_point_x, center_point_y):
        distance = 0
        for projected_point in projected_points:
          distance += ((projected_point[0] - center_point_x)**2 + (projected_point[1] - center_point_y)**2)**-4
        
        # add penalty if optimized point is too far from cog
        distance_to_cog = (center_point_x - center_point[0])**2+(center_point_y - center_point[1])**2
        
        distance += distance_to_cog*1e-8   
        return distance
      
      # compute the rotation angle when iterating over all connection vectors between center and border point
      def rotation_angle(center_point_x, center_point_y):
        total_angle = 0
        last_vector = None
        for projected_point in projected_points:
          vector = [-center_point_x + projected_point[0], -center_point_y + projected_point[1]]
          #print("projected_point: {}, center_point: ({},{}), vector: {}, last_vector: {}".format(projected_point, center_point_x, center_point_y, vector, last_vector))
          if last_vector is not None:
            denominator = np.sqrt(vector[0]**2 + vector[1]**2) * np.sqrt(last_vector[0]**2 + last_vector[1]**2)
            value = -(vector[0] * last_vector[1] - vector[1] * last_vector[0])/denominator
            angle = np.arcsin(value)
            #print("value: ", value, ", angle: ", angle*180./np.pi)
            total_angle += angle
            
          last_vector = list(vector)
          
        projected_point = projected_points[0,:]
        vector = [-center_point_x + projected_point[0], -center_point_y + projected_point[1]]
        #print("first projected_point: {}, center_point: ({},{}), vector: {}, last_vector: {}".format(projected_point, center_point_x, center_point_y, vector, last_vector))
        denominator = np.sqrt(vector[0]**2 + vector[1]**2) * np.sqrt(last_vector[0]**2 + last_vector[1]**2)
        value = -(vector[0] * last_vector[1] - vector[1] * last_vector[0])/denominator
        angle = np.arcsin(value)
        #print("angle: ", angle*180./np.pi)
        total_angle += angle
            
        return total_angle
        
      #a = rotation_angle(center_point[0], center_point[1])
      #print("test a=",a*180./np.pi)
      
      import casadi   # sudo pip install casadi

      # Symbols/expressions
      x = casadi.MX.sym('x')
      y = casadi.MX.sym('y')
      f = squared_distance_to_all_points(x,y)

      nlp = {}                 # NLP declaration
      nlp['x']= casadi.vertcat(x,y) # decision vars
      nlp['f'] = f             # objective
      #nlp['g'] = rotation_angle(x,y)             # constraints
      previous_center_point = [center_point[0], center_point[1]]
      initial_values = [center_point[0], center_point[1]]

      # Create solver instance
      F = casadi.nlpsol('F','ipopt',nlp);

      # Solve the problem using a guess
      #result = F(x0=initial_values, ubg=370./180.*np.pi, lbg=350./180.*np.pi)
      result = F(x0=initial_values)
      center_point[0] = result['x'][0]
      center_point[1] = result['x'][1]
      #print("previous_center_point: ", previous_center_point, ", optimized center point: ", center_point)
    
      a = rotation_angle(center_point[0], center_point[1])
      #print("resulting rotation_angle: ",a*180./np.pi)
      distance_to_cog = np.sqrt((previous_center_point[0] - center_point[0])**2+(previous_center_point[1] - center_point[1])**2)
      #print("resulting distance to cog: ", distance_to_cog)
        
    # add center point as new point
    projected_points = np.concatenate([projected_points, np.array([[center_point[0],center_point[1]]])],axis=0)
    
    # transform projected points back to 3D points
    points = []
    for projected_point in projected_points:
      points.append(np.array([projected_point[0], projected_point[1], z_value]))
    
    # update n_points
    n_points = len(points)
    points = np.reshape(points, (-1,3))
    
    center_point_index = n_points-1
    
    # create list with point indices for each triangle
    point_indices_list = []
    for i in range(len(points)-2):
      point_indices_list.append([center_point_index, i, i+1])
    point_indices_list.append([center_point_index, n_points-2, 0])
    
    
    #print("point_indices_list:",point_indices_list)
    #print("points:",points,points.shape)
    triangle_list = []
    for point_index_triple in point_indices_list:
      triangle_list.append(points[point_index_triple])
    
    #triangle_list = points[point_indices_list]  ## doesn't work sometimes
    
  print("  number of projected points: ",len(projected_points),", number of initial triangles: ", len(point_indices_list))
    
  # solve Laplace equation on deformed mesh
  n_dofs = n_points
  global_stiffness = np.zeros((n_dofs, n_dofs))
  
  # assemble global stiffness matrix
  for point_indices in point_indices_list:
    
    tri = points[point_indices]
    
    # extract points as column vectors
    point0 = np.reshape(tri[0],(-1,1))
    point1 = np.reshape(tri[1],(-1,1))
    point2 = np.reshape(tri[2],(-1,1))
    
    if debugging_stl_output:
      out_triangulation_world_space.append([tri[0], tri[1], tri[2]])
    
    # compute Jacobian of x(xi) mapping
    # xi-space: triangle (0,0) (1,0) (0,1),     (parameter space)
    # x-space: triangle point0 point1 point2, x^{1}, x^{2}, x^{3}    (given by initial triangulation)
    # u-space: triangle u^{1} u^{2} u^{3}   (to be computed)
    # barycentric coordinates
    # x(xi), xi(u)
    # barycentric coordinates: 
    #   x(xi) = (1-xi1-xi2)*x^{1} + xi1*x^{2} + xi2*x^{3}
    #   u(xi) = (1-xi1-xi2)*u^{1} + xi1*u^{2} + xi2*u^{3}
    #   u(xi) = u  =>  (1-xi1-xi2)*u^{1}_1 + xi1*u^{2}_1 + xi2*u^{3}_1 = u
    #              =>  (1-xi1-xi2)*u^{1}_2 + xi1*u^{2}_2 + xi2*u^{3}_2 = v
    #
    #                  xi1*(u^{2}_1 - u^{1}_1)  +  xi2*(u^{3}_1 - u^{1}_1)  =  u - u^{1}_1
    #                  xi1*(u^{2}_2 - u^{1}_2)  +  xi2*(u^{3}_2 - u^{1}_2)  =  u - u^{1}_2
    #   xi(u) = u_{,xi}^-1 (u - u^{1}) with u_{,xi} = (u^{2}-u^{1}  u^{3} - u^{1}), and xi_{,u} = u_{,xi}^-1
    #
    #   x(u) = x(xi(u))
    #   x_{,u} = x_{,xi} * xi_{,u}
    #
    #   M_x(u) = x_{,u}^T x_{,u}
    #   M_x(xi) = x_{,xi}^T x_{,xi}
    
    dxdxi = np.concatenate([point1-point0, point2-point0],axis=1)
    metric_dxdxi = np.transpose(dxdxi).dot(dxdxi)   # dxdxi^T * dxdxi
    metric_dxdxi_inv = np.linalg.inv(metric_dxdxi)
    integration_factor = np.sqrt(np.linalg.det(metric_dxdxi))
    
    def phi(xi1,xi2,dof_no):
      # dof no:
      #    2
      #    |\
      #    |_\
      #   0   1
      if dof_no == 0:
        return (1. - xi1) * (1. - xi2)
      
      elif dof_no == 1:
        return xi1 * (1. - xi2)
        
      elif dof_no == 2:
        return (1. - xi1) * xi2
        
    def grad_phi(xi1,xi2,dof_no):
      if dof_no == 0:
        return np.array([[-(1. - xi2)], [-(1. - xi1)]])
      
      elif dof_no == 1:
        return np.array([[(1. - xi2)], [-xi1]])
        
      elif dof_no == 2:
        return np.array([[-xi2], [(1. - xi1)]])
        
    # compute element stiffness matrix
    def integrand(xi1,xi2,i,j):
      result = np.transpose(grad_phi(xi1,xi2,i)).dot(metric_dxdxi_inv).dot(grad_phi(xi1,xi2,j)) * integration_factor
      return result
        
    element_stiffness = np.empty((3,3))
    for i in range(3):
      for j in range(3):
        result = scipy.integrate.dblquad(integrand, 0.0, 1.0, lambda xi1: 0.0, lambda xi2: 1.0-xi2, args=(i,j))
        element_stiffness[i][j] = result[0]
      
    # assemble global stiffness matrix
    for i in range(3):
      for j in range(3):
        dof_i = point_indices[i]
        dof_j = point_indices[j]
        global_stiffness[dof_i][dof_j] += element_stiffness[i][j]
        
  
  # assemble rhs vector of Dirchlet BC
  dirichlet_u = np.zeros((n_dofs,1))
  dirichlet_v = np.zeros((n_dofs,1))
  rhs_u = np.zeros((n_dofs,1))
  rhs_v = np.zeros((n_dofs,1))
  
  if debug:
    print("n_dofs=",n_dofs,",n_original_points=",n_original_points)
  
  # loop over boundary points, `original_points` are the points of the ring surface, `points` is a superset containing additional points created by the triangulation
  for original_point_no,original_point in enumerate(original_points):
    
    # get the position in reference space
    if parametric_space_shape == 0 or parametric_space_shape == 3 or parametric_space_shape == 4:  # unit circle
      phi = float(original_point_no) / n_original_points * 2 * np.pi
      if triangulation_type == 1:
        phi = original_point_phi_value[original_point_no]
      u_reference = np.cos(phi)
      v_reference = np.sin(phi)
      
    elif parametric_space_shape == 1:  # unit square
      if original_point_no < n_original_points/4:   # bottom
        u_reference = float(original_point_no)/(n_original_points/4)
        v_reference = 0.
      elif original_point_no < 2*n_original_points/4:   # right
        u_reference = 1.
        v_reference = float(original_point_no-n_original_points/4)/(n_original_points/4)
      elif original_point_no < 3*n_original_points/4:   # top
        u_reference = 1.0 - float(original_point_no-2*n_original_points/4)/(n_original_points/4)
        v_reference = 1.
      else:   # left
        u_reference = 0.0
        v_reference = 1.0 - float(original_point_no-3*n_original_points/4)/(n_original_points/4)
    
    elif parametric_space_shape == 2:  # unit square with adjusted grid position
      if original_point_no < n_original_points/4:   # bottom        
        alpha = float(original_point_no)/(n_original_points/4)
        phi = -np.pi/4 + np.pi/2.*alpha
        a = 0.5    # distance between center of square and line of points
        
        u_reference = 0.5 + np.tan(phi)*a
        v_reference = 0.
        
      elif original_point_no < 2*n_original_points/4:   # right        
        alpha = float(original_point_no-n_original_points/4)/(n_original_points/4)
        phi = -np.pi/4 + np.pi/2.*alpha
        a = 0.5    # distance between center of square and line of points
        
        u_reference = 1.
        v_reference = 0.5 + np.tan(phi)*a
        
      elif original_point_no < 3*n_original_points/4:   # top        
        alpha = 1.0 - float(original_point_no-2*n_original_points/4)/(n_original_points/4)
        phi = -np.pi/4 + np.pi/2.*alpha
        a = 0.5    # distance between center of square and line of points
        
        u_reference = 0.5 + np.tan(phi)*a
        v_reference = 1.
        
      else:   # left
        alpha = 1.0 - float(original_point_no-3*n_original_points/4)/(n_original_points/4)
        phi = -np.pi/4 + np.pi/2.*alpha
        a = 0.5    # distance between center of square and line of points
        
        u_reference = 0.0
        v_reference = 0.5 + np.tan(phi)*a
    
    # find the corresponding point in points
    dof_no = None
    for no,point in enumerate(points):
      if abs(point[0]-original_point[0]) < 1e-12 and abs(point[1]-original_point[1]) < 1e-12:
        dof_no = no
        break
        
    if dof_no is not None:
      dirichlet_u[dof_no] = u_reference
      dirichlet_v[dof_no] = v_reference

  # There are 2 types of techniques to ensure the Dirichlet BC: one explict and one implicit type.
  # Both are implemented and both yield the same result

  # strong form Dirichlet BC
  if False:
    # loop over dirichlet bc dof index
    for i in range(n_original_points):
      
      # loop over all dofs
      for j in range(n_dofs):
        if i == j:
          rhs_u[i] = dirichlet_u[i]
          rhs_v[i] = dirichlet_v[i]
        else:
          rhs_u[j] -= global_stiffness[j][i] * dirichlet_u[i]
          rhs_v[j] -= global_stiffness[j][i] * dirichlet_v[i]


    # loop over dirichlet bc dof index
    for i in range(n_original_points):
      rhs_u[i] = dirichlet_u[i]
      rhs_v[i] = dirichlet_v[i]
    
    # modify stiffness matrix for Dirichlet boundary conditions
    global_stiffness[0:n_original_points,:] = 0.0
    global_stiffness[:,0:n_original_points] = 0.0
    for i in range(n_original_points):
      global_stiffness[i,i] = 1.0

    if debug:
      print("global stiffness matrix:")
      print(global_stiffness)

      print("")
      print("rhs vector (u):")
      print(rhs_u)

    # solve liner system
    u = np.linalg.solve(global_stiffness, rhs_u)
    v = np.linalg.solve(global_stiffness, rhs_v)
    
  # weak form Dirichlet BC
  if True:

    # Dirichlet indices are 0:n_original_points
    # non-Dirichlet indices are n_original_points:

    global_stiffness_non_dirichlet = global_stiffness[n_original_points:,n_original_points:]

    # compute rhs vector    
    global_stiffness_dirichlet = global_stiffness[n_original_points:,0:n_original_points]
    rhs_non_dirichlet_u = -global_stiffness_dirichlet.dot(dirichlet_u[0:n_original_points])
    rhs_non_dirichlet_v = -global_stiffness_dirichlet.dot(dirichlet_v[0:n_original_points])

    # solve liner system
    non_dirichlet_u = np.linalg.solve(global_stiffness_non_dirichlet, rhs_non_dirichlet_u)
    non_dirichlet_v = np.linalg.solve(global_stiffness_non_dirichlet, rhs_non_dirichlet_v)
    
    if debug:
      print("rhs_non_dirichlet_u:",rhs_non_dirichlet_u)
      print("global_stiffness_non_dirichlet:",global_stiffness_non_dirichlet)
      print("non_dirichlet_u:",non_dirichlet_u)
      print("non_dirichlet_v:",non_dirichlet_v)
    
    u = np.concatenate([dirichlet_u[0:n_original_points], non_dirichlet_u])
    v = np.concatenate([dirichlet_v[0:n_original_points], non_dirichlet_v])
    
    if debug:
      print("u:",u)
    
  if debug:
    # output solution
    for dof_no in range(n_dofs):
      
      # find the corresponding point in original_points
      point = points[dof_no]
      original_point_no = None
      for no,original_point in enumerate(original_points):
        if abs(point[0]-original_point[0]) < 1e-12 and abs(point[1]-original_point[1]) < 1e-12:
          original_point_no = no
          break
        
      print("dof {}, original point no: {}, dirichlet: ({},{}), solution: ({},{}), rhs: ({},{})".\
        format(dof_no, original_point_no, dirichlet_u[dof_no], dirichlet_v[dof_no], u[dof_no], v[dof_no], rhs_u[dof_no], rhs_v[dof_no]))
      
  # store the triangles in parametric space
  triangles_parametric_space = []
  
  # loop over triangles indices
  for point_indices in point_indices_list:
    
    # loop over the dofs of the current triangle and store (u,v) values of current triangle
    triangle_parametric_space = []
    for dof_no in point_indices:
      point_parametric_space = np.array([u[dof_no], v[dof_no]])
      
      triangle_parametric_space.append(point_parametric_space)
    
    # store triangle to list of triangles
    triangles_parametric_space.append(triangle_parametric_space)
    
  if debugging_stl_output:
    # create parametric space triangles for debugging output, move near ring and scale 10x
    x_offset = center_point[0] + extent_x*1.5
    y_offset = center_point[1]
    scale = 10.0
    for triangle_parametric_space in triangles_parametric_space:
      
      out_triangle_parametric_space = []
      for p in triangle_parametric_space:
        out_triangle_parametric_space.append(np.array([p[0]*scale+x_offset, p[1]*scale+y_offset, z_value]))
        
      out_triangulation_parametric_space.append(out_triangle_parametric_space)
      
  # now the mapping x -> u,v is computed (from parameter space to world space)
  # create new grid points on the ring that form a uniform mesh in parametric space

  n_grid_points = n_grid_points_x*n_grid_points_y
  grid_points_world_space = np.empty((n_grid_points,3))
  grid_points_parametric_space = np.empty((n_grid_points,2))
  grid_points_parametric_space_modified = np.empty((n_grid_points,2))
  
  # loop over grid points in parametric space
  for (j,y) in enumerate(np.linspace(0.0,1.0,n_grid_points_y)):
    phi = float(y) * (n_grid_points_y-1.0) / n_grid_points_y  * 2.*np.pi
    for (i,x) in enumerate(np.linspace(0.0,1.0,n_grid_points_x)):
  
      if parametric_space_shape == 0:  # unit circle
        r = x
        x = r*np.cos(phi)
        y = r*np.sin(phi)
      
      elif parametric_space_shape == 2:  # square with adjusted grid points
        # get segment
        if (n_grid_points_y%2 == 1 and j != int(n_grid_points_y/2.)) or n_grid_points_y%2 == 0:
          if j < n_grid_points_y/2.:   # bottom or side
            if i >= j and i <= n_grid_points_x-1-j:   # bottom quarter /\
          
              # get layer
              a = ((n_grid_points_y-1)/2.-j)*1.0/(n_grid_points_y-1)    # distance between center of square and line of points
              alpha = float(i - j) / (n_grid_points_x-1-2*j)
              phi = -np.pi/4. + np.pi/2.*alpha
              
              x = 0.5 + np.tan(phi)*a
          else:   # top or side
            if i >= n_grid_points_x-1-j and i <= j:   # top
          
              # get layer
              a = (j-(n_grid_points_y-1)/2.)*1.0/(n_grid_points_y-1)    # distance between center of square and line of points
              alpha = float(i - (n_grid_points_y-1-j)) / (2*j-n_grid_points_y+1)
              phi = -np.pi/4 + np.pi/2.*alpha
              
              x = 0.5 + np.tan(phi)*a
            
        if (n_grid_points_x%2 == 1 and i != int(n_grid_points_x/2.)) or n_grid_points_y%2 == 0:
          if i < n_grid_points_x/2.:   # left
            if j >= i and j <= n_grid_points_y-1-i:   # left
          
              # get layer
              a = ((n_grid_points_x-1)/2.-i)*1.0/(n_grid_points_x-1)    # distance between center of square and line of points
              alpha = float(j - i) / (n_grid_points_y-1-2*i)
              phi = -np.pi/4 + np.pi/2.*alpha
              
              y = 0.5 + np.tan(phi)*a
          else:   # right
            if j >= n_grid_points_y-i and j <= i:   # right
          
              # get layer
              a = (i-(n_grid_points_x-1)/2.)*1.0/(n_grid_points_x-1)    # distance between center of square and line of points
              alpha = float(j - (n_grid_points_y-1-i)) / (2*i-n_grid_points_x+1)
              phi = -np.pi/4 + np.pi/2.*alpha
              
              y = 0.5 + np.tan(phi)*a
        
      elif parametric_space_shape == 3 or parametric_space_shape == 4:    # unit circle with adjusted grid points
        
        # get segment
        if (n_grid_points_y%2 == 1 and j != int(n_grid_points_y/2.)) or n_grid_points_y%2 == 0:
          if j < n_grid_points_y/2.:   # bottom or side
            if i >= j and i <= n_grid_points_x-1-j:   # bottom
          
              # get layer
              a = ((n_grid_points_y-1)/2.-j)*2.0/(n_grid_points_y-1)    # fraction (in [0,1]) of distance between center of square and line of points
              alpha = float(i - j) / (n_grid_points_x-1-2*j)
              phi = -np.pi/4. + np.pi/2.*alpha
              
              x = np.sin(phi)*a
              y = (-1./np.sqrt(2.) + (-np.cos(phi) + 1./np.sqrt(2.))*a)*a
          
          else:   # top or side
            if i >= n_grid_points_x-1-j and i <= j:   # top
          
              # get layer
              a = (j-(n_grid_points_y-1)/2.)*2.0/(n_grid_points_y-1)    # distance between center of square and line of points
              alpha = float(i - (n_grid_points_y-1-j)) / (2*j-n_grid_points_y+1)
              phi = -np.pi/4 + np.pi/2.*alpha
              
              x = np.sin(phi)*a
              y = (1./np.sqrt(2.) + (np.cos(phi) - 1./np.sqrt(2.))*a)*a
          
        if (n_grid_points_x%2 == 1 and i != int(n_grid_points_x/2.)) or n_grid_points_y%2 == 0:
          if i < n_grid_points_x/2.:   # left
            if j >= i and j <= n_grid_points_y-1-i:   # left
          
              # get layer
              a = ((n_grid_points_x-1)/2.-i)*2.0/(n_grid_points_x-1)    # distance between center of square and line of points
              alpha = float(j - i) / (n_grid_points_y-1-2*i)
              phi = -np.pi/4 + np.pi/2.*alpha
              
              y = np.sin(phi)*a
              x = (-1./np.sqrt(2.) + (-np.cos(phi) + 1./np.sqrt(2.))*a)*a
          
          else:   # right
            if j >= n_grid_points_y-i and j <= i:   # right
          
              # get layer
              a = (i-(n_grid_points_x-1)/2.)*2.0/(n_grid_points_x-1)    # distance between center of square and line of points
              alpha = float(j - (n_grid_points_y-1-i)) / (2*i-n_grid_points_x+1)
              phi = -np.pi/4 + np.pi/2.*alpha
              
              y = np.sin(phi)*a
              x = (1./np.sqrt(2.) + (np.cos(phi) - 1./np.sqrt(2.))*a)*a
           
        # rotate by pi*3/4 
        phi = np.arctan2(y, x) + np.pi/4 + np.pi/2
        r = np.sqrt(x*x + y*y)
        x = np.cos(phi)*r
        y = np.sin(phi)*r
        
        # modified
        phi_modified = get_modified_phi(phi)
        x_modified = np.cos(phi_modified)*r
        y_modified = np.sin(phi_modified)*r
        
        if n_grid_points_x%2 == 1 and i == int(n_grid_points_x/2) and j == int(n_grid_points_y/2):   # center point
          x = 0.
          y = 0.
          x_modified = 0.
          y_modified = 0.
          
      # transform to world space
      # if there were additional border points from the triangulation that also got Dirichlet BC, and modifiy_phi is set to true, compute world points differently
      if modify_phi:
        point_world_space = transform_to_world_space(x_modified,y_modified,triangles_parametric_space,triangle_list,parametric_space_shape)
      else:
        point_world_space = transform_to_world_space(x,y,triangles_parametric_space,triangle_list,parametric_space_shape)
    
      if debug and False:
        print("modify_phi: {} [{},{}] -> {} (triangles_parametric_space={}, triangle_list={}, parametric_space_shape={})".
          format(modify_phi, x,y, point_world_space, triangles_parametric_space, triangle_list, parametric_space_shape))
    
      if point_world_space is None:
        grid_points_world_space[j*n_grid_points_x+i] = np.array([0.0,0.0,0.0])
    
      if point_world_space is not None:
        grid_points_world_space[j*n_grid_points_x+i] = point_world_space
        grid_points_parametric_space[j*n_grid_points_x+i] = np.array([x,y])
        
        if modify_phi:
          grid_points_parametric_space_modified[j*n_grid_points_x+i] = np.array([x_modified,y_modified])
        else:
          grid_points_parametric_space_modified[j*n_grid_points_x+i] = grid_points_parametric_space[j*n_grid_points_x+i]
    
  if parametric_space_shape == 4:
    # move grid points such that mean distance between points in world space gets optimal
    
    # objective function
    def objective(inner_grid_points_parametric):
          
      inner_grid_points_parametric = np.reshape(inner_grid_points_parametric, (-1,2))
      #print("shape of inner_grid_points_parametric: ", inner_grid_points_parametric.shape
      
      grid_points_world = np.zeros((n_grid_points_y*n_grid_points_x,3))
        
      # transform grid points from parametric space to world space
      for j in range(0,n_grid_points_y):
        for i in range(0,n_grid_points_x):
          if i > 0 and i < n_grid_points_x-1 and j > 0 and j < n_grid_points_y-1:
            
            p = inner_grid_points_parametric[(j-1)*(n_grid_points_x-2)+(i-1),:]
            point_world = transform_to_world_space(p[0],p[1],triangles_parametric_space,triangle_list,parametric_space_shape)
          
            if point_world is None:
              grid_points_world[j*n_grid_points_x+i,:] = np.array([0,0,0])
          
            if point_world is not None:
              grid_points_world[j*n_grid_points_x+i,:] = np.array([point_world[0],point_world[1],z_value])
            
          else:
            grid_points_world[j*n_grid_points_x+i,:] = grid_points_world_space[j*n_grid_points_x+i,:]
          
      #print("grid_points_world:",grid_points_world
      
      distances_current_loop,relative_distances_current_loop = compute_mean_distances(grid_points_world, n_grid_points_x, n_grid_points_y)
      #print("objective, std: {}".format(np.std(relative_distances_current_loop))
      return np.std(relative_distances_current_loop)
    
    
    # set initial values
    initial_values = np.zeros(((n_grid_points_y-2)*(n_grid_points_x-2),2))
    
    for j in range(0,n_grid_points_y-2):
      for i in range(0,n_grid_points_x-2):
        initial_values[j*(n_grid_points_x-2)+i,:] = grid_points_parametric_space[(j+1)*n_grid_points_x+(i+1),:]
    initial_values = np.reshape(initial_values, (-1,1))
      
    #print("initial_values: ",initial_values
      
    result = scipy.optimize.minimize(objective, initial_values, method='Nelder-Mead', options={"maxiter":1e4, "disp":True})
    print(result["message"])
    resulting_parametric_points = result["x"]
    
    print("final objective: {}".format(objective(resulting_parametric_points)))
    
    for j in range(0,n_grid_points_y):
      for i in range(0,n_grid_points_x):
        if i > 0 and i < n_grid_points_x-1 and j > 0 and j < n_grid_points_y-1:
          grid_points_parametric_space[j*n_grid_points_x+i,0] = resulting_parametric_points[2*((j-1)*(n_grid_points_x-2)+(i-1))+0]
          grid_points_parametric_space[j*n_grid_points_x+i,1] = resulting_parametric_points[2*((j-1)*(n_grid_points_x-2)+(i-1))+1]
      
        x = grid_points_parametric_space[j*n_grid_points_x+i,0]
        y = grid_points_parametric_space[j*n_grid_points_x+i,1]
      
        # transform to world space
        point_world_space = transform_to_world_space(x,y,triangles_parametric_space,triangle_list,parametric_space_shape)
      
        #print("point_world_space:",point_world_space
      
        if point_world_space is None:
          #print("point_world_space is None"
          grid_points_world_space[j*n_grid_points_x+i] = np.array([0.0,0.0,0.0])
      
        if point_world_space is not None:
          grid_points_world_space[j*n_grid_points_x+i] = point_world_space
          grid_points_parametric_space[j*n_grid_points_x+i] = np.array([x,y])
    

  #distances_current_loop,relative_distances_current_loop = compute_mean_distances(grid_points_world_space)
  #print("transformed, std: {}".format(np.std(relative_distances_current_loop))
  
  # improve grid_points_world_space
  grid_points_world_space_improved = grid_points_world_space
  
  if improve_mesh:   # optimize points, fix self-intersections and Laplacian smoothing
    import copy
    import random
    random.seed(0)
    print("  improving mesh (fixing self-intersections and Laplacian smoothing)")
    
    grid_points_world_space_improved = copy.deepcopy(grid_points_world_space)
    
  
    def ccw(p0,p1,p2):
      """ if triangle p0,p1,p2 is counterclockwise, source: https://algs4.cs.princeton.edu/91primitives/ """
      #print((p1[0]-p0[0])*(p2[1]-p0[1]) - (p2[0]-p0[0])*(p1[1]-p0[1]))
      return ((p1[0]-p0[0])*(p2[1]-p0[1]) - (p2[0]-p0[0])*(p1[1]-p0[1])) > 0.0
    
    def ccw_value(p0,p1,p2):
      """ if triangle p0,p1,p2 is counterclockwise, source: https://algs4.cs.princeton.edu/91primitives/ """
      #print((p1[0]-p0[0])*(p2[1]-p0[1]) - (p2[0]-p0[0])*(p1[1]-p0[1]))
      return ((p1[0]-p0[0])*(p2[1]-p0[1]) - (p2[0]-p0[0])*(p1[1]-p0[1]))
    
    def lines_intersect(ap,aq,bp,bq):
      """ check if line (ap,aq) intersects line (bp,bq) """
      if (ccw_value(ap, aq, bp) * ccw_value(ap, aq, bq) > 0):
        return False
      if (ccw_value(bp, bq, ap) * ccw_value(bp, bq, aq) > 0):
        return False
      return True

    def is_self_intersecting(p0,p1,p2,p3):
      """ if the quadrilateral is self-intersecting """
      # p2 p3
      # p0 p1
      
      return lines_intersect(p0,p1,p2,p3) or lines_intersect(p0,p2,p1,p3)
    
    #def is_self_intersecting(p0,p1,p2,p3):
    #  """ if the quadrilateral is self-intersecting """
    #  # p2 p3
    #  # p0 p1
    #  return (not ccw(p0,p1,p3) or not ccw(p1,p3,p2) or not ccw(p3,p2,p0) or not ccw(p2,p0,p1))
    
    def quadrilateral_contains_point(p0,p1,p2,p3,p):
      return triangle_contains_point([p0,p1,p2],p)[0] or triangle_contains_point([p2,p1,p3],p)[0]
      
    def is_properly_oriented(p0,p1,p2,p3):
      """ if quadrilateral is oriented counter-clockwise, this also excludes is_self_intersecting """
      # p2 p3
      # p0 p1
      n_ccw = (1 if ccw(p0,p1,p3) else 0) + (1 if ccw(p1,p3,p2) else 0) + (1 if ccw(p3,p2,p0) else 0) + (1 if ccw(p2,p0,p1) else 0)
      
      return n_ccw >= 3
      
    def does_overlap(p0,p1p2,p3p4):
      """ check if the triangles [p2 p0 p1] and [p4 p0 p3] overlap """
      [p1,p2] = p1p2
      [p3,p4] = p3p4
      a14 = oriented_angle(-p0+p1,-p0+p4) 
      a13 = oriented_angle(-p0+p1,-p0+p3) 
      a12 = oriented_angle(-p0+p1,-p0+p2) 
      a23 = oriented_angle(-p0+p2,-p0+p3) 
      a24 = oriented_angle(-p0+p2,-p0+p4) 
      a31 = oriented_angle(-p0+p3,-p0+p1) 
      a34 = oriented_angle(-p0+p3,-p0+p4)
      a32 = oriented_angle(-p0+p3,-p0+p2)
      
      
      if (a13 > 0 and a12 > 0 and a13 < a12 and abs(a13) > 1e-5 and abs(a23) > 1e-5)\
        or (a13 < 0 and a12 < 0 and a13 > a12 and abs(a13) > 1e-5 and abs(a23) > 1e-5):
        return True
      if (a14 > 0 and a12 > 0 and a14 < a12 and abs(a14) > 1e-5 and abs(a24) > 1e-5)\
        or (a14 < 0 and a12 < 0 and a14 > a12 and abs(a14) > 1e-5 and abs(a24) > 1e-5):
        return True
           
      if (a31 > 0 and a34 > 0 and a31 < a34 and abs(a31) > 1e-5 and abs(a14) > 1e-5)\
        or (a31 < 0 and a34 < 0 and a31 > a34 and abs(a31) > 1e-5 and abs(a14) > 1e-5):
        return True
      if (a32 > 0 and a34 > 0 and a32 < a34 and abs(a32) > 1e-5 and abs(a24) > 1e-5)\
        or (a32 < 0 and a34 < 0 and a32 > a34 and abs(a32) > 1e-5 and abs(a24) > 1e-5):
        return True
        
      return False
    
    def oriented_angle(v1,v2):
      """ compute the angle between v1 and v2, with proper sign """
      phi = np.arccos(v1.dot(v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
      orientation = v1[0]*v2[1] - v1[1]*v2[0]
      if abs(orientation) < 1e-12:
        if abs(v2[0]) > abs(v2[1]):
          if v1[0]/v2[0] > 0:
            return 0
          else:
            return np.pi
        else:
          if v1[1]/v2[1] > 0:
            return 0
          else:
            return np.pi
      elif orientation > 0:
        return phi
      else:
        return -phi
      
      return angle
    
    # compute score of elements, a low score means a good shaped element
    def get_element_score(p0,p1,p2,p3):
      
      debug = False
      
      # variance of side lengths, favours quadrilaterals with same side lengths
      p01 = p1-p0
      p12 = p3-p1
      p23 = p2-p3
      p30 = p0-p2
      characteristic_length = (extent_x+extent_y)/2.
      l01 = np.linalg.norm(p01) / characteristic_length
      l12 = np.linalg.norm(p12) / characteristic_length
      l23 = np.linalg.norm(p23) / characteristic_length
      l30 = np.linalg.norm(p30) / characteristic_length
      
      #variance_side_lengths = np.var([l01, l12, l23, l30]) * 5
      
      cutoff = 10  # lower cutoff means penalize short side lengths stronger, which leads to longer side lengths
      variance_side_lengths = max(0, 1./l01-cutoff) + max(0, 1./l12-cutoff) + max(0, 1./l23-cutoff) + max(0, 1./l30-cutoff)   # penalize small side lengths
      variance_side_lengths *= 5e-2
      
      # angles
      a0 = np.arctan2(np.linalg.norm(np.cross(p01, -p30)), np.dot(p01, -p30))
      a1 = np.arctan2(np.linalg.norm(np.cross(p12, -p01)), np.dot(p12, -p01))
      a2 = np.arctan2(np.linalg.norm(np.cross(p23, -p12)), np.dot(p23, -p12))
      a3 = np.arctan2(np.linalg.norm(np.cross(p30, -p23)), np.dot(p30, -p23))
      
      
      #angles_score1 = np.linalg.norm([a0-np.pi/4., a1-np.pi/4., a2-np.pi/4., a3-np.pi/4.])   # penalize derivation from 90 degrees
      
      angles_score2 = max(0.,1./(a0**1)-1) + max(0.,1./(a1**1)-1) + max(0.,1./(a2**1)-1) + max(0.,1./(a3**1)-1)   # penalize small angles
      
      #angles_score3 = (max(0, -a0) + max(0, -a1) + max(0, -a2) + max(0, -a3))*10000  # penalize negative angles
      
      # massively penalize self-intersecting quadrilaterals
      angles_score3 = 0 if is_properly_oriented(p0,p1,p2,p3) else 10000
      
      #print("ccw: ",ccw(p0,p1,p3), ccw(p1,p3,p2), ccw(p3,p2,p0), ccw(p2,p0,p1))
      
      #angles_score4 = np.var([a0,a1,a2,a3])*1e-1
      angles_score4 = 0
      
      #print(angles_score2, angles_score4)
      variance_angles_lengths = (angles_score2 + angles_score3 + angles_score4)
      
      score = variance_angles_lengths + variance_side_lengths
      
      if debug:
        print("side_lengths: {} {} {} {} score: {}, angles: {} {} {}, score: {}, total: {}".format(1./l01, 1./l12, 1./l23, 1./l30, variance_side_lengths, angles_score2, angles_score3, angles_score4, variance_angles_lengths, score))
      return score
      
    # compute score for all elements
    if False:  # not needed      
      n_elements = (n_grid_points_x-1) * (n_grid_points_y-1)
      element_score = np.zeros(n_elements)
      for j in range(0,n_grid_points_y-1):
        for i in range(0,n_grid_points_x-1):
          # p2 p3
          # p0 p1
            
          p0 = grid_points_world_space_improved[j*n_grid_points_x+i]
          p1 = grid_points_world_space_improved[j*n_grid_points_x+(i+1)%n_grid_points_x]
          p2 = grid_points_world_space_improved[(j+1)%n_grid_points_y*n_grid_points_x+i]
          p3 = grid_points_world_space_improved[(j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x]
          
          score = get_element_score(p0,p1,p2,p3)
          element_score[j*(n_grid_points_x-1)+i] = score
          
         
    # output grid
    output_pre_fix = False
    output_fix = False
    output_post_fix = False
    if output_pre_fix:
      patches_world_improved = []
      
      # loop over grid points in parametric space
      for (jj,y) in enumerate(np.linspace(0.0,1.0,n_grid_points_y)):
        phi = float(y) * (n_grid_points_y-1.0) / n_grid_points_y  * 2.*np.pi
        for (ii,x) in enumerate(np.linspace(0.0,1.0,n_grid_points_x)):
    
          if parametric_space_shape == 1 or parametric_space_shape == 2 or parametric_space_shape == 3 or parametric_space_shape == 4:  # unit square  
            if ii == n_grid_points_x-1 or jj == n_grid_points_x-1:
              continue
          if parametric_space_shape == 0:
            if ii == n_grid_points_x-1:
              continue
            
          p0_improved = grid_points_world_space_improved[jj*n_grid_points_x+ii]
          p1_improved = grid_points_world_space_improved[jj*n_grid_points_x+(ii+1)%n_grid_points_x]
          p2_improved = grid_points_world_space_improved[(jj+1)%n_grid_points_y*n_grid_points_x+ii]
          p3_improved = grid_points_world_space_improved[(jj+1)%n_grid_points_y*n_grid_points_x+(ii+1)%n_grid_points_x]
          
          quadrilateral = np.zeros((4,2))
          quadrilateral[0] = p0_improved[0:2]
          quadrilateral[1] = p1_improved[0:2]
          quadrilateral[2] = p3_improved[0:2]
          quadrilateral[3] = p2_improved[0:2]
          
          min_x = min(min_x, min(quadrilateral[:,0]))
          min_y = min(min_y, min(quadrilateral[:,1]))
          max_x = max(max_x, max(quadrilateral[:,0]))
          max_y = max(max_y, max(quadrilateral[:,1]))
          
          polygon = patches.Polygon(quadrilateral, True)
          patches_world_improved.append(polygon)
          
      # world space, improved
      fig, ax = plt.subplots(figsize=(20,20))
        
      xw_improved = np.reshape(grid_points_world_space_improved[:,0], (-1))
      yw_improved = np.reshape(grid_points_world_space_improved[:,1], (-1))

      patch_collection = collections.PatchCollection(patches_world_improved,edgecolors="k",facecolors="gray",alpha=0.5)
      ax.add_collection(patch_collection)
      ax.plot(xw_improved, yw_improved, "ok")
      ax.set_xlim(min_x,max_x)
      ax.set_ylim(min_y,max_y)
      plt.axis('equal')
      
      plt.savefig("out/loop_{:03}_p{}_world_mesh_pre_fix.png".format(loop_no, os.getpid()));
      if show_plot:
        plt.show()
      plt.close()
        
    factor = (extent_x*extent_y)/2 * 5e-3
    
    # try to resolve self-intersecting quadrilaterals
    # repeatedly iterate over all elements until no more fixes could be achieved
    while True:
      changed_a_point = False
      are_all_elements_properly_oriented = True
      
      # loop over all elements
      for i in range(0,n_grid_points_x-1):
        for j in range(0,n_grid_points_y-1):
          # p2 p3
          # p0 p1
      
          p0 = grid_points_world_space_improved[j*n_grid_points_x+i]
          p1 = grid_points_world_space_improved[j*n_grid_points_x+(i+1)%n_grid_points_x]
          p2 = grid_points_world_space_improved[(j+1)%n_grid_points_y*n_grid_points_x+i]
          p3 = grid_points_world_space_improved[(j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x]
    
          if not is_properly_oriented(p0,p1,p2,p3):
            are_all_elements_properly_oriented = False
            if output_fix:
              print("  self intersection: p0=np.array([{},{},{}]); p1=np.array([{},{},{}]); p2=np.array([{},{},{}]); p3=np.array([{},{},{}]) # {}+{}+{}+{}<3".format(p0[0],p0[1],p0[2],p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],p3[0],p3[1],p3[2],(1 if ccw(p0,p1,p3) else 0),(1 if ccw(p1,p3,p2) else 0),(1 if ccw(p3,p2,p0) else 0),(1 if ccw(p2,p0,p1) else 0)))
            else:
              print("  self intersection in element ({},{})".format(i,j))
            
            indices = [(i,j),(i+1,j),(i,j+1),(i+1,j+1)]
            
            # loop over 4 points of element
            for k in range(4):
              (ii,jj) = indices[k]
              
              # do not consider border points, they cannot be changed
              if ii <= 0 or jj <= 0 or ii >= n_grid_points_x-1 or jj >= n_grid_points_y-1:
                continue
              
              if output_fix:
                print("({},{})/({},{})".format(ii,jj,n_grid_points_x,n_grid_points_y))
              
              # the point p will be moved, p0-p7 are the neighbouring points
              p = grid_points_world_space_improved[jj*n_grid_points_x+ii]
              p_old = np.array(p)
              # p6 p5 p4
              # p7 p  p3
              # p0 p1 p2
              p0 = grid_points_world_space_improved[(jj-1)*n_grid_points_x+(ii-1)]
              p1 = grid_points_world_space_improved[(jj-1)*n_grid_points_x+ii]
              p2 = grid_points_world_space_improved[(jj-1)*n_grid_points_x+(ii+1)]
              p3 = grid_points_world_space_improved[jj*n_grid_points_x+(ii+1)]
              p4 = grid_points_world_space_improved[(jj+1)*n_grid_points_x+(ii+1)]
              p5 = grid_points_world_space_improved[(jj+1)*n_grid_points_x+ii]
              p6 = grid_points_world_space_improved[(jj+1)*n_grid_points_x+(ii-1)]
              p7 = grid_points_world_space_improved[jj*n_grid_points_x+(ii-1)]
            
              if output_fix:
                plt.figure(figsize=(20,20))
                              
                plt.plot(p[0],p[1],'go')
                plt.plot(p0[0],p0[1],'ro')
                plt.plot(p1[0],p1[1],'ro')
                plt.plot(p2[0],p2[1],'ro')
                plt.plot(p3[0],p3[1],'ro')
                plt.plot(p4[0],p4[1],'ro')
                plt.plot(p5[0],p5[1],'ro')
                plt.plot(p6[0],p6[1],'ro')
                plt.plot(p7[0],p7[1],'ro')

                plt.plot([p0[0],p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0],p0[0]], [p0[1],p1[1],p2[1],p3[1],p4[1],p5[1],p6[1],p7[1],p0[1]], 'k-')
                plt.plot([p1[0],p[0]], [p1[1],p[1]], 'k-')
                plt.plot([p3[0],p[0]], [p3[1],p[1]], 'k-')
                plt.plot([p5[0],p[0]], [p5[1],p[1]], 'k-')
                plt.plot([p7[0],p[0]], [p7[1],p[1]], 'k-')
                  
                s = ""
                if does_overlap(p,[p1,p3],[p7,p5]):
                  s += "e"
                if not is_properly_oriented(p0,p1,p7,p):
                  s += "f"
                if not is_properly_oriented(p1,p2,p,p3):
                  s += "g"
                if not is_properly_oriented(p7,p,p6,p5):
                  s += "h"
                if not is_properly_oriented(p,p3,p5,p4):
                  s += "i"
                
                plt.savefig("out/{}_{}_areference_out_{}_{}.png".format(i,j,k,s))
                        
              n_tries = 0
              size_factor = factor
              p_changed = np.array(p)
              
              def orientation_score(p0,p1,p2,p3):
                """ how well the quad is oriented, 0=worst, 3=properly """
                # p2 p3
                # p0 p1
                n_ccw = (1 if ccw(p0,p1,p3) else 0) + (1 if ccw(p1,p3,p2) else 0) + (1 if ccw(p3,p2,p0) else 0) + (1 if ccw(p2,p0,p1) else 0)
                return min(n_ccw,3)
              
              initial_score = orientation_score(p0,p1,p7,p_changed) + orientation_score(p1,p2,p_changed,p3) \
                + orientation_score(p7,p_changed,p6,p5) + orientation_score(p_changed,p3,p5,p4)
              
              # while the score is not yet better (and maximum of 200 tries), deflect point p
              while (\
                (orientation_score(p0,p1,p7,p_changed) + orientation_score(p1,p2,p_changed,p3) \
                + orientation_score(p7,p_changed,p6,p5) + orientation_score(p_changed,p3,p5,p4)) <= initial_score \
                and n_tries < 200):
                
                # pseudo-randomly deflect point p
                p_changed = p + np.array([(random.random()-0.5)*size_factor, (random.random()-0.5)*size_factor, 0])
                size_factor *= 1.05    # 1.05**200 = 17292
                
                if output_fix:
                  plt.figure()
                                
                  plt.plot(p[0],p[1],'ko')
                  plt.plot(p_changed[0],p_changed[1],'ro')
                  plt.plot(p0[0],p0[1],'bo')
                  plt.plot(p1[0],p1[1],'bo')
                  plt.plot(p2[0],p2[1],'bo')
                  plt.plot(p3[0],p3[1],'bo')
                  plt.plot(p4[0],p4[1],'bo')
                  plt.plot(p5[0],p5[1],'bo')
                  plt.plot(p6[0],p6[1],'bo')
                  plt.plot(p7[0],p7[1],'bo')

                  plt.plot([p0[0],p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0],p0[0]], [p0[1],p1[1],p2[1],p3[1],p4[1],p5[1],p6[1],p7[1],p0[1]], 'k-')
                  plt.plot([p1[0],p_changed[0]], [p1[1],p_changed[1]], 'k-')
                  plt.plot([p3[0],p_changed[0]], [p3[1],p_changed[1]], 'k-')
                  plt.plot([p5[0],p_changed[0]], [p5[1],p_changed[1]], 'k-')
                  plt.plot([p7[0],p_changed[0]], [p7[1],p_changed[1]], 'k-')
                  
                  s = ""
                  #if does_overlap(p,[p1,p3],[p7,p5]):
                  #  s += "e"
                  if not is_properly_oriented(p0,p1,p7,p_changed):
                    s += "f"
                  if not is_properly_oriented(p1,p2,p_changed,p3):
                    s += "g"
                  if not is_properly_oriented(p7,p_changed,p6,p5):
                    s += "h"
                  if not is_properly_oriented(p_changed,p3,p5,p4):
                    s += "i"
                  
                  plt.savefig("out/{}_{}_out_{}_{}_{}_{}.png".format(i,j,k,n_tries,size_factor,s))
                  plt.close()
                  
                n_tries += 1
              
              if n_tries < 200:
                if is_properly_oriented(p0,p1,p2,p3):
                  print("  \033[0;32mSuccessfully resolved self-intersection ({},{}) after {} iterations\033[0m".format(i,j,n_tries))
                else:
                  current_score = (orientation_score(p0,p1,p7,p_changed) + orientation_score(p1,p2,p_changed,p3) \
                    + orientation_score(p7,p_changed,p6,p5) + orientation_score(p_changed,p3,p5,p4))
                  print("  improvement regarding self-intersection ({},{}) after {} iterations, score: {} -> {}".format(i,j,n_tries, initial_score, current_score))
                grid_points_world_space_improved[jj*n_grid_points_x+ii] = p_changed
                changed_a_point = True
                
                if debugging_stl_output:
                  plt.figure()
                                
                  plt.plot(p[0],p[1],'ko')
                  plt.plot(p_changed[0],p_changed[1],'ro')
                  plt.plot(p0[0],p0[1],'bo')
                  plt.plot(p1[0],p1[1],'bo')
                  plt.plot(p2[0],p2[1],'bo')
                  plt.plot(p3[0],p3[1],'bo')
                  plt.plot(p4[0],p4[1],'bo')
                  plt.plot(p5[0],p5[1],'bo')
                  plt.plot(p6[0],p6[1],'bo')
                  plt.plot(p7[0],p7[1],'bo')

                  plt.plot([p0[0],p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0],p0[0]], [p0[1],p1[1],p2[1],p3[1],p4[1],p5[1],p6[1],p7[1],p0[1]], 'k-')
                  plt.plot([p1[0],p_changed[0]], [p1[1],p_changed[1]], 'k-')
                  plt.plot([p3[0],p_changed[0]], [p3[1],p_changed[1]], 'k-')
                  plt.plot([p5[0],p_changed[0]], [p5[1],p_changed[1]], 'k-')
                  plt.plot([p7[0],p_changed[0]], [p7[1],p_changed[1]], 'k-')
                  
                  score = (orientation_score(p0,p1,p7,p_changed) + orientation_score(p1,p2,p_changed,p3) \
                    + orientation_score(p7,p_changed,p6,p5) + orientation_score(p_changed,p3,p5,p4))
                  
                  plt.savefig("out/{}_{}_out_{}_{}_{}_{}.png".format(i,j,k,n_tries,size_factor,score))
                  plt.close()
                  
                
                break
              else:
                print("  \033[0;31mself-intersection was not resolved after {} iterations\033[0m".format(n_tries))
            
            if output_fix:
              p0 = grid_points_world_space_improved[j*n_grid_points_x+i]
              p1 = grid_points_world_space_improved[j*n_grid_points_x+(i+1)%n_grid_points_x]
              p2 = grid_points_world_space_improved[(j+1)%n_grid_points_y*n_grid_points_x+i]
              p3 = grid_points_world_space_improved[(j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x]
              print("p0=np.array({}); p1=np.array({}); p2=np.array({}); p3=np.array({}); ".format(p0,p1,p2,p3))
              
            
        #print("n_grid_points_x: {}, n_grid_points_y: {}, size: {}".format(n_grid_points_x, n_grid_points_y, len(grid_points_world_space_improved)))
      
      # if there was a resolved self_intersection, restart iterations, otherwise we are done
      if not changed_a_point:
        if are_all_elements_properly_oriented:
          print("  all elements are properly oriented")
        break
      
    # output grid
    if output_post_fix:
      patches_world_improved = []
      
      # loop over grid points in parametric space
      for (jj,y) in enumerate(np.linspace(0.0,1.0,n_grid_points_y)):
        phi = float(y) * (n_grid_points_y-1.0) / n_grid_points_y  * 2.*np.pi
        for (ii,x) in enumerate(np.linspace(0.0,1.0,n_grid_points_x)):
    
          if parametric_space_shape == 1 or parametric_space_shape == 2 or parametric_space_shape == 3 or parametric_space_shape == 4:  # unit square  
            if ii == n_grid_points_x-1 or jj == n_grid_points_x-1:
              continue
          if parametric_space_shape == 0:
            if ii == n_grid_points_x-1:
              continue
              
          p0_improved = grid_points_world_space_improved[jj*n_grid_points_x+ii]
          p1_improved = grid_points_world_space_improved[jj*n_grid_points_x+(ii+1)%n_grid_points_x]
          p2_improved = grid_points_world_space_improved[(jj+1)%n_grid_points_y*n_grid_points_x+ii]
          p3_improved = grid_points_world_space_improved[(jj+1)%n_grid_points_y*n_grid_points_x+(ii+1)%n_grid_points_x]
          
          quadrilateral = np.zeros((4,2))
          quadrilateral[0] = p0_improved[0:2]
          quadrilateral[1] = p1_improved[0:2]
          quadrilateral[2] = p3_improved[0:2]
          quadrilateral[3] = p2_improved[0:2]
          
          min_x = min(min_x, min(quadrilateral[:,0]))
          min_y = min(min_y, min(quadrilateral[:,1]))
          max_x = max(max_x, max(quadrilateral[:,0]))
          max_y = max(max_y, max(quadrilateral[:,1]))
          
          polygon = patches.Polygon(quadrilateral, True)
          patches_world_improved.append(polygon)
          
      # world space, improved
      fig, ax = plt.subplots(figsize=(20,20))
        
      xw_improved = np.reshape(grid_points_world_space_improved[:,0], (-1))
      yw_improved = np.reshape(grid_points_world_space_improved[:,1], (-1))

      patch_collection = collections.PatchCollection(patches_world_improved,edgecolors="k",facecolors="gray",alpha=0.5)
      ax.add_collection(patch_collection)
      ax.plot(xw_improved, yw_improved, "ok")
      ax.set_xlim(min_x,max_x)
      ax.set_ylim(min_y,max_y)
      plt.axis('equal')
      
      plt.savefig("out/loop_{:03}_p{}_world_mesh_post_fix.png".format(loop_no, os.getpid()));
      if show_plot:
        plt.show()
      plt.close()
    
    # improve point locations by Laplacian smoothing
    random.seed(1)
    for k in range(20):
      
      for i in range(1,n_grid_points_x-1):
        for j in range(1,n_grid_points_y-1):
                    
          p = grid_points_world_space_improved[j*n_grid_points_x+i]
          p_old = np.array(p)
          # p6 p5 p4
          # p7 p  p3
          # p0 p1 p2
          p0 = grid_points_world_space_improved[(j-1)*n_grid_points_x+(i-1)]
          p1 = grid_points_world_space_improved[(j-1)*n_grid_points_x+i]
          p2 = grid_points_world_space_improved[(j-1)*n_grid_points_x+(i+1)]
          p3 = grid_points_world_space_improved[j*n_grid_points_x+(i+1)]
          p4 = grid_points_world_space_improved[(j+1)*n_grid_points_x+(i+1)]
          p5 = grid_points_world_space_improved[(j+1)*n_grid_points_x+i]
          p6 = grid_points_world_space_improved[(j+1)*n_grid_points_x+(i-1)]
          p7 = grid_points_world_space_improved[j*n_grid_points_x+(i-1)]
        
          p_changed = 0.125 * (p0+p1+p2+p3+p4+p5+p6+p7)
          
          grid_points_world_space_improved[j*n_grid_points_x+i] = p_changed
          
          # output grid
          output = False
          if output:
            # print output in objective function
            objective(p_changed,True)
            
            patches_world_improved = []
            
            # loop over grid points in parametric space
            for (jj,y) in enumerate(np.linspace(0.0,1.0,n_grid_points_y)):
              phi = float(y) * (n_grid_points_y-1.0) / n_grid_points_y  * 2.*np.pi
              for (ii,x) in enumerate(np.linspace(0.0,1.0,n_grid_points_x)):
          
                if parametric_space_shape == 1 or parametric_space_shape == 2 or parametric_space_shape == 3 or parametric_space_shape == 4:  # unit square  
                  if ii == n_grid_points_x-1 or jj == n_grid_points_x-1:
                    continue
                if parametric_space_shape == 0:
                  if ii == n_grid_points_x-1:
                    continue
                    
                p0_improved = grid_points_world_space_improved[jj*n_grid_points_x+ii]
                p1_improved = grid_points_world_space_improved[jj*n_grid_points_x+(ii+1)%n_grid_points_x]
                p2_improved = grid_points_world_space_improved[(jj+1)%n_grid_points_y*n_grid_points_x+ii]
                p3_improved = grid_points_world_space_improved[(jj+1)%n_grid_points_y*n_grid_points_x+(ii+1)%n_grid_points_x]
                
                quadrilateral = np.zeros((4,2))
                quadrilateral[0] = p0_improved[0:2]
                quadrilateral[1] = p1_improved[0:2]
                quadrilateral[2] = p3_improved[0:2]
                quadrilateral[3] = p2_improved[0:2]
                
                min_x = min(min_x, min(quadrilateral[:,0]))
                min_y = min(min_y, min(quadrilateral[:,1]))
                max_x = max(max_x, max(quadrilateral[:,0]))
                max_y = max(max_y, max(quadrilateral[:,1]))
                
                polygon = patches.Polygon(quadrilateral, True)
                patches_world_improved.append(polygon)
                
            # world space, improved
            fig, ax = plt.subplots(figsize=(20,20))
              
            xw_improved = np.reshape(grid_points_world_space_improved[:,0], (-1))
            yw_improved = np.reshape(grid_points_world_space_improved[:,1], (-1))

            patch_collection = collections.PatchCollection(patches_world_improved,edgecolors="k",facecolors="gray",alpha=0.5)
            ax.add_collection(patch_collection)
            ax.plot(xw_improved, yw_improved, "ok")
            ax.plot(p_old[0],p_old[1], 'rx')
            ax.plot(p_changed[0],p_changed[1], 'g+')
            ax.set_xlim(min_x,max_x)
            ax.set_ylim(min_y,max_y)
            plt.axis('equal')
            
            plt.savefig("out/{}{}{}_loop_{:03}_p{}_world_mesh_improved_k{}_{}_{}.png".format(k,i,j,loop_no, os.getpid(), k, old_score, new_score));
            if show_plot:
              plt.show()
            plt.close()
    
  if debugging_stl_output:
    # create triangles of new grid points mesh
    grid_point_indices_world_space = []
    
    patches_parametric = []
    patches_world = []
    patches_world_improved = []
    patches_parametric_modified = []
    parametric_points = []
    parametric_points_modified = []
    min_x = 100000
    min_y = 100000
    max_x = -100000
    max_y = -100000
    factor = 1.0
    
    # loop over grid points in parametric space
    for (j,y) in enumerate(np.linspace(0.0,1.0,n_grid_points_y)):
      phi = float(y) * (n_grid_points_y-1.0) / n_grid_points_y  * 2.*np.pi
      for (i,x) in enumerate(np.linspace(0.0,1.0,n_grid_points_x)):
    
        [x,y] = grid_points_parametric_space[j*n_grid_points_x+i]
    
        # for debugging create markers at grid points in parametric space
        point = np.array(np.array([x*scale+x_offset, y*scale+y_offset, z_value]))
        factor *= 1.04
        if factor >= 3:
          factor = 3
        size = 0.1*factor
        diag0 = np.array([-size,-size,-size])
        diag1 = np.array([size,-size,-size])
        diag2 = np.array([-size,size,-size])
        diag3 = np.array([size,size,-size])
        diag4 = np.array([-size,-size,size])
        diag5 = np.array([size,-size,size])
        diag6 = np.array([-size,size,size])
        diag7 = np.array([size,size,size])
        markers_grid_points_parametric_space += [
          [point+diag0,point+diag3,point+diag1],[point+diag0,point+diag2,point+diag3],  # bottom
          [point+diag4,point+diag5,point+diag7],[point+diag4,point+diag7,point+diag6],  # top
          [point+diag0,point+diag1,point+diag5],[point+diag0,point+diag5,point+diag4],  # front
          [point+diag2,point+diag7,point+diag3],[point+diag2,point+diag6,point+diag7],  # back
          [point+diag2,point+diag0,point+diag4],[point+diag2,point+diag4,point+diag6],  # left
          [point+diag1,point+diag3,point+diag7],[point+diag1,point+diag7,point+diag5]  # right
        ]

        # for debugging create markers at grid points in world space
        point = grid_points_world_space[j*n_grid_points_x+i]
        size = 0.1*factor
        diag0 = np.array([-size,-size,-size])
        diag1 = np.array([size,-size,-size])
        diag2 = np.array([-size,size,-size])
        diag3 = np.array([size,size,-size])
        diag4 = np.array([-size,-size,size])
        diag5 = np.array([size,-size,size])
        diag6 = np.array([-size,size,size])
        diag7 = np.array([size,size,size])
        markers_grid_points_world_space += [
          [point+diag0,point+diag3,point+diag1],[point+diag0,point+diag2,point+diag3],  # bottom
          [point+diag4,point+diag5,point+diag7],[point+diag4,point+diag7,point+diag6],  # top
          [point+diag0,point+diag1,point+diag5],[point+diag0,point+diag5,point+diag4],  # front
          [point+diag2,point+diag7,point+diag3],[point+diag2,point+diag6,point+diag7],  # back
          [point+diag2,point+diag0,point+diag4],[point+diag2,point+diag4,point+diag6],  # left
          [point+diag1,point+diag3,point+diag7],[point+diag1,point+diag7,point+diag5]  # right
        ]

        if parametric_space_shape == 1 or parametric_space_shape == 2 or parametric_space_shape == 3 or parametric_space_shape == 4:  # unit square  
          if i == n_grid_points_x-1 or j == n_grid_points_x-1:
            continue
        if parametric_space_shape == 0:
          if i == n_grid_points_x-1:
            continue
        # p2 p3
        # p0 p1
          
        p0 = grid_points_world_space[j*n_grid_points_x+i]
        p1 = grid_points_world_space[j*n_grid_points_x+(i+1)%n_grid_points_x]
        p2 = grid_points_world_space[(j+1)%n_grid_points_y*n_grid_points_x+i]
        p3 = grid_points_world_space[(j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x]
        
        grid_point_indices_world_space.append([j*n_grid_points_x+i, j*n_grid_points_x+(i+1)%n_grid_points_x, (j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x])
        grid_point_indices_world_space.append([j*n_grid_points_x+i, (j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x, (j+1)%n_grid_points_y*n_grid_points_x+i])
        
        grid_triangles_world_space.append([p0,p1,p3])
        grid_triangles_world_space.append([p0,p3,p2])
        
        p0_improved = grid_points_world_space_improved[j*n_grid_points_x+i]
        p1_improved = grid_points_world_space_improved[j*n_grid_points_x+(i+1)%n_grid_points_x]
        p2_improved = grid_points_world_space_improved[(j+1)%n_grid_points_y*n_grid_points_x+i]
        p3_improved = grid_points_world_space_improved[(j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x]
        
        quadrilateral = np.zeros((4,2))
        quadrilateral[0] = p0[0:2]
        quadrilateral[1] = p1[0:2]
        quadrilateral[2] = p3[0:2]
        quadrilateral[3] = p2[0:2]
        
        min_x = min(min_x, min(quadrilateral[:,0]))
        min_y = min(min_y, min(quadrilateral[:,1]))
        max_x = max(max_x, max(quadrilateral[:,0]))
        max_y = max(max_y, max(quadrilateral[:,1]))
        #print("world: ",quadrilateral
        polygon = patches.Polygon(quadrilateral, True)
        patches_world.append(polygon)
        
        quadrilateral = np.zeros((4,2))
        quadrilateral[0] = p0_improved[0:2]
        quadrilateral[1] = p1_improved[0:2]
        quadrilateral[2] = p3_improved[0:2]
        quadrilateral[3] = p2_improved[0:2]
        polygon = patches.Polygon(quadrilateral, True)
        patches_world_improved.append(polygon)
        
        offset = np.array([x_offset, y_offset])
        p0 = np.concatenate([grid_points_parametric_space[j*n_grid_points_x+i]*scale+offset,np.array([z_value])])
        p1 = np.concatenate([grid_points_parametric_space[j*n_grid_points_x+(i+1)%n_grid_points_x]*scale+offset,np.array([z_value])])
        p2 = np.concatenate([grid_points_parametric_space[(j+1)%n_grid_points_y*n_grid_points_x+i]*scale+offset,np.array([z_value])])
        p3 = np.concatenate([grid_points_parametric_space[(j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x]*scale+offset,np.array([z_value])])
        
        grid_triangles_parametric_space.append([p0,p1,p3])
        grid_triangles_parametric_space.append([p0,p3,p2])
        
        quadrilateral = np.zeros((4,2))
        quadrilateral[0] = grid_points_parametric_space[j*n_grid_points_x+i]
        quadrilateral[1] = grid_points_parametric_space[j*n_grid_points_x+(i+1)%n_grid_points_x]
        quadrilateral[2] = grid_points_parametric_space[(j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x]
        quadrilateral[3] = grid_points_parametric_space[(j+1)%n_grid_points_y*n_grid_points_x+i]
        parametric_points.append(quadrilateral[0])
        parametric_points.append(quadrilateral[1])
        parametric_points.append(quadrilateral[2])
        parametric_points.append(quadrilateral[3])
        #print("parametric: ",quadrilateral
        polygon = patches.Polygon(quadrilateral, True)
        patches_parametric.append(polygon)
        
        quadrilateral = np.zeros((4,2))
        quadrilateral[0] = grid_points_parametric_space_modified[j*n_grid_points_x+i]
        quadrilateral[1] = grid_points_parametric_space_modified[j*n_grid_points_x+(i+1)%n_grid_points_x]
        quadrilateral[2] = grid_points_parametric_space_modified[(j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x]
        quadrilateral[3] = grid_points_parametric_space_modified[(j+1)%n_grid_points_y*n_grid_points_x+i]
        parametric_points_modified.append(quadrilateral[0])
        parametric_points_modified.append(quadrilateral[1])
        parametric_points_modified.append(quadrilateral[2])
        parametric_points_modified.append(quadrilateral[3])
        #print("parametric: ",quadrilateral
        polygon = patches.Polygon(quadrilateral, True)
        patches_parametric_modified.append(polygon)
        
  t_stop = timeit.default_timer()
  duration = t_stop - t_start

  if debugging_stl_output:
    # plot laplace solutions
    x = np.reshape(points[:,0], (-1))
    y = np.reshape(points[:,1], (-1))
    
    x0 = x[0:n_points_per_face]
    x1 = x[n_points_per_face:2*n_points_per_face]
    x2 = x[2*n_points_per_face:3*n_points_per_face]
    x3 = x[3*n_points_per_face:4*n_points_per_face]
    x4 = x[4*n_points_per_face:4*n_points_per_face+n_additional_points_on_ring]
    x5 = x[4*n_points_per_face+n_additional_points_on_ring:]
    
    y0 = y[0:n_points_per_face]
    y1 = y[n_points_per_face:2*n_points_per_face]
    y2 = y[2*n_points_per_face:3*n_points_per_face]
    y3 = y[3*n_points_per_face:4*n_points_per_face]
    y4 = y[4*n_points_per_face:4*n_points_per_face+n_additional_points_on_ring]
    y5 = y[4*n_points_per_face+n_additional_points_on_ring:]
    
    u_list = np.reshape(u, (-1))
    v_list = np.reshape(v, (-1))
    
    u0 = u_list[0:n_points_per_face]
    u1 = u_list[n_points_per_face:2*n_points_per_face]
    u2 = u_list[2*n_points_per_face:3*n_points_per_face]
    u3 = u_list[3*n_points_per_face:4*n_points_per_face]
    u4 = u_list[4*n_points_per_face:4*n_points_per_face+n_additional_points_on_ring]
    u5 = u_list[4*n_points_per_face+n_additional_points_on_ring:]
    
    v0 = v_list[0:n_points_per_face]
    v1 = v_list[n_points_per_face:2*n_points_per_face]
    v2 = v_list[2*n_points_per_face:3*n_points_per_face]
    v3 = v_list[3*n_points_per_face:4*n_points_per_face]
    v4 = v_list[4*n_points_per_face:4*n_points_per_face+n_additional_points_on_ring]
    v5 = v_list[4*n_points_per_face+n_additional_points_on_ring:]
    
    xw = np.reshape(grid_points_world_space[:,0], (-1))
    yw = np.reshape(grid_points_world_space[:,1], (-1))

    xw_improved = np.reshape(grid_points_world_space_improved[:,0], (-1))
    yw_improved = np.reshape(grid_points_world_space_improved[:,1], (-1))

    f, ax = plt.subplots(2,3,figsize=(30,20))
    
    # u
    ax[0,0].tricontourf(x,y,u_list, 100) # 100 contour levels
    ax[0,0].triplot(x,y,point_indices_list,color='k')
    ax[0,0].plot(x0,y0, 'ro')
    ax[0,0].plot(x1,y1, 'yo')
    ax[0,0].plot(x2,y2, 'go')
    ax[0,0].plot(x3,y3, 'bo')
    ax[0,0].plot(x4,y4, 'md')
    ax[0,0].plot(x5,y5, 'ko')
    ax[0,0].set_title('u, triangulation in world space')
    ax[0,0].set_aspect('equal')
    
    # v
    ax[0,1].tricontourf(x,y,v_list, 100) # 100 contour levels
    ax[0,1].triplot(x,y,point_indices_list,color='k')
    ax[0,1].plot(x0,y0, 'ro')
    ax[0,1].plot(x1,y1, 'yo')
    ax[0,1].plot(x2,y2, 'go')
    ax[0,1].plot(x3,y3, 'bo')
    ax[0,1].plot(x4,y4, 'md')
    ax[0,1].plot(x5,y5, 'ko')
    ax[0,1].set_title('v, triangulation in world space')
    ax[0,1].set_aspect('equal')
    
    # parametric space
    if determine_additional_points_on_ring:
      p = collections.PatchCollection(patches_parametric_modified,edgecolors="k",facecolors="white")
      ax[0,2].add_collection(p)
      ax[0,2].plot([p[0] for p in parametric_points_modified],[p[1] for p in parametric_points_modified], 'ko')
      if parametric_space_shape == 1 or parametric_space_shape == 2:
        ax[0,2].set_xlim(-0.1,1.1)
        ax[0,2].set_ylim(-0.1,1.1)
      else:
        ax[0,2].set_xlim(-1.1,1.1)
        ax[0,2].set_ylim(-1.1,1.1)
      plt.axis('equal')
        
      for j in range(n_grid_points_y):
        for i in range(n_grid_points_x):
          p = grid_points_parametric_space_modified[j*n_grid_points_x+i]
          if j == 0:
            ax[0,2].plot(p[0], p[1], 'ro')
          elif j == n_grid_points_y-1:
            ax[0,2].plot(p[0], p[1], 'go')
          elif i == 0: 
            ax[0,2].plot(p[0], p[1], 'bo')
          elif i == n_grid_points_x-1: 
            ax[0,2].plot(p[0], p[1], 'yo')
          else:
            ax[0,2].plot(p[0], p[1], 'ko')
    else:
      p = collections.PatchCollection(patches_parametric,edgecolors="k",facecolors="white")
      ax[0,2].add_collection(p)
      ax[0,2].plot([p[0] for p in parametric_points],[p[1] for p in parametric_points], 'ko')
      
      for j in range(n_grid_points_y):
        for i in range(n_grid_points_x):
          p = grid_points_parametric_space[j*n_grid_points_x+i]
          if j == 0:
            ax[0,2].plot(p[0], p[1], 'ro')
          elif j == n_grid_points_y-1:
            ax[0,2].plot(p[0], p[1], 'go')
          elif i == 0: 
            ax[0,2].plot(p[0], p[1], 'bo')
          elif i == n_grid_points_x-1: 
            ax[0,2].plot(p[0], p[1], 'yo')
          else:
            ax[0,2].plot(p[0], p[1], 'ko')
        
    ax[0,2].set_title('quadrangulation in parametric space')
    if parametric_space_shape == 1 or parametric_space_shape == 2:
      ax[0,2].set_xlim(-0.1,1.1)
      ax[0,2].set_ylim(-0.1,1.1)
    else:
      ax[0,2].set_xlim(-1.1,1.1)
      ax[0,2].set_ylim(-1.1,1.1)
    
    # parametric space
    ax[1,0].triplot(u_list,v_list,point_indices_list,color='k')
    ax[1,0].plot(u0,v0, 'ro')
    ax[1,0].plot(u1,v1, 'yo')
    ax[1,0].plot(u2,v2, 'go')
    ax[1,0].plot(u3,v3, 'bo')
    ax[1,0].plot(u4,v4, 'md')
    ax[1,0].plot(u5,v5, 'ko')
    ax[1,0].set_title('triangulation in parametric space')
    ax[1,0].set_aspect('equal')
    
    # world space grid
    ax[1,1].triplot(xw,yw,grid_point_indices_world_space,color='k')
    
    for j in range(n_grid_points_y):
      for i in range(n_grid_points_x):
        p = grid_points_world_space[j*n_grid_points_x+i]
        if j == 0:
          ax[1,1].plot(p[0], p[1], 'ro')
        elif j == n_grid_points_y-1:
          ax[1,1].plot(p[0], p[1], 'go')
        elif i == 0: 
          ax[1,1].plot(p[0], p[1], 'bo')
        elif i == n_grid_points_x-1: 
          ax[1,1].plot(p[0], p[1], 'yo')
        else:
          ax[1,1].plot(p[0], p[1], 'ko')
          
    ax[1,1].set_title('new grid in world space')
    ax[1,1].set_aspect('equal')
    
    # world space grid
    ax[1,2].triplot(xw_improved,yw_improved,grid_point_indices_world_space,color='k')
    for j in range(n_grid_points_y):
      for i in range(n_grid_points_x):
        p = grid_points_world_space_improved[j*n_grid_points_x+i]
        if j == 0:
          ax[1,2].plot(p[0], p[1], 'ro')
        elif j == n_grid_points_y-1:
          ax[1,2].plot(p[0], p[1], 'go')
        elif i == 0: 
          ax[1,2].plot(p[0], p[1], 'bo')
        elif i == n_grid_points_x-1: 
          ax[1,2].plot(p[0], p[1], 'yo')
        else:
          ax[1,2].plot(p[0], p[1], 'ko')
          
    ax[1,2].set_title('after improving')
    ax[1,2].set_aspect('equal')
    
    filename = "out/loop_{:03}_p{}_harmonic_map.png".format(loop_no, os.getpid())
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
      print("Create directory \"{}\".".format(dirname))
      os.makedirs(dirname)
    print("Save \"{}\".".format(filename))
    
    matplotlib.rcParams.update({
      'font.size': 24,
      'lines.markersize': 10
    })
      
    plt.savefig(filename)
    if show_plot:
      plt.show()
    plt.close()
    
    # plot quadrilaterals
    if False:
      # parametric space
      fig, ax = plt.subplots()
        
      p = collections.PatchCollection(patches_parametric,edgecolors="k",facecolors="gray",alpha=0.5)
      ax.add_collection(p)
      ax.plot([p[0] for p in parametric_points],[p[1] for p in parametric_points], 'ko')
      ax.set_xlim(-1.1,1.1)
      ax.set_ylim(-1.1,1.1)
      plt.axis('equal')
      
      plt.savefig("out/loop_{:03}_p{}_parametric_mesh.png".format(loop_no, os.getpid()));
      if show_plot:
        plt.show()
      plt.close()
    
    # parametric space with modified phi
    if determine_additional_points_on_ring:
      fig, ax = plt.subplots(figsize=(10,10))
    
      p = collections.PatchCollection(patches_parametric_modified,edgecolors="k",facecolors="white")
      ax.add_collection(p)
      ax.plot([p[0] for p in parametric_points_modified],[p[1] for p in parametric_points_modified], 'ko')
      if parametric_space_shape == 1 or parametric_space_shape == 2:
        ax.set_xlim(-0.1,1.1)
        ax.set_ylim(-0.1,1.1)
      else:
        ax.set_xlim(-1.1,1.1)
        ax.set_ylim(-1.1,1.1)
      plt.axis('equal')
        
      for j in range(n_grid_points_y):
        for i in range(n_grid_points_x):
          p = grid_points_parametric_space_modified[j*n_grid_points_x+i]
          if j == 0:
            ax.plot(p[0], p[1], 'ro')
          elif j == n_grid_points_y-1:
            ax.plot(p[0], p[1], 'go')
          elif i == 0: 
            ax.plot(p[0], p[1], 'bo')
          elif i == n_grid_points_x-1: 
            ax.plot(p[0], p[1], 'yo')
          else:
            ax.plot(p[0], p[1], 'ko')
      plt.savefig("out/loop_{:03}_p{}_modified_parametric_mesh.png".format(loop_no, os.getpid()));
      if show_plot:
        plt.show()
      plt.close()
      
    # world space
    fig, ax = plt.subplots(figsize=(20,20))
      
    patch_collection_world = collections.PatchCollection(patches_world,edgecolors="k",facecolors="gray",alpha=0.5)
    ax.add_collection(patch_collection_world)
    #ax.plot(xw,yw, 'ko',markersize=10)
    ax.plot(xw,yw, 'ko')
    ax.set_xlim(min_x,max_x)
    ax.set_ylim(min_y,max_y)
    plt.axis('equal')
    
    plt.savefig("out/loop_{:03}_p{}_world_mesh.png".format(loop_no, os.getpid()));
    if show_plot:
      plt.show()
    plt.close()
    
    # world space, improved
    fig, ax = plt.subplots(figsize=(20,20))
      
    patch_collection_improved = collections.PatchCollection(patches_world_improved,edgecolors="k",facecolors="gray",alpha=0.5)
    ax.add_collection(patch_collection_improved)
    #ax.plot(xw,yw, 'ko',markersize=10)
    ax.plot(xw,yw, 'yx')
    ax.plot(xw_improved,yw_improved, 'ko')
    ax.set_xlim(min_x,max_x)
    ax.set_ylim(min_y,max_y)
    plt.axis('equal')
    
    plt.savefig("out/loop_{:03}_p{}_world_mesh_improved.png".format(loop_no, os.getpid()));
    if show_plot:
      plt.show()
    plt.close()
       
  return grid_points_world_space_improved, duration
 
def create_3d_mesh(loop_grid_points, n_grid_points_x, n_grid_points_y, debugging_stl_output, out_3d_mesh_triangles):
  """
  Create a 3D mesh by connecting the 2d meshes.
  :param loop_grid_points: List of list of the grid points for each loop
  :param debugging_stl_output: if list should be filled with STL triangles that can be output to a STL mesh for debugging
  :param out_3d_mesh_triangles: output list of triangles 
  """
  
  # create 3D mesh from grid points on slices
  node_positions = []
  nodes = []
  linear_elements = []   # 8 node indices per element
  quadratic_elements = []   # 27 node indices per element

  # check if also quadratic elements are possible from the number of grid points
  create_quadratic_elements = True
  n_loops = len(loop_grid_points)
  if n_grid_points_x <= 2 or n_grid_points_x%2 == 0 or n_grid_points_y <= 2 or n_grid_points_y%2 == 0 or n_loops%2 == 0 or n_loops <= 2:
    create_quadratic_elements = False
    print("Grid per ring is {} x {}, {} rings: with this it is not possible to create quadratic elements (numbers must be odd)".format(n_grid_points_x, n_grid_points_y, n_loops))

  # fill list of nodes
  for grid_point_list in loop_grid_points:
    for node_position in list(grid_point_list):
      node_positions.append(node_position.tolist())
    nodes += list(grid_point_list)
    
  # fill list of elements
  n_grid_points_per_loop = n_grid_points_x * n_grid_points_y
  for loop_no in range(n_loops-1):
    
    grid_point_list_bottom = loop_grid_points[loop_no]
    grid_point_list_top = loop_grid_points[loop_no+1]
    
    # linear elements
    for j in range(n_grid_points_y-1):
      for i in range(n_grid_points_x-1):
        
        i0 = loop_no*n_grid_points_per_loop + j*n_grid_points_x+i
        i1 = loop_no*n_grid_points_per_loop + j*n_grid_points_x+i+1
        i2 = loop_no*n_grid_points_per_loop + (j+1)*n_grid_points_x+i
        i3 = loop_no*n_grid_points_per_loop + (j+1)*n_grid_points_x+i+1
        i4 = (loop_no+1)*n_grid_points_per_loop + j*n_grid_points_x+i
        i5 = (loop_no+1)*n_grid_points_per_loop + j*n_grid_points_x+i+1
        i6 = (loop_no+1)*n_grid_points_per_loop + (j+1)*n_grid_points_x+i
        i7 = (loop_no+1)*n_grid_points_per_loop + (j+1)*n_grid_points_x+i+1
        
        linear_elements.append([i0,i1,i2,i3,i4,i5,i6,i7])
        
  for loop_no in range(n_loops-2):
    # quadratic elements
    if create_quadratic_elements:
      for j in range(n_grid_points_y-2):
        for i in range(n_grid_points_x-2):
          
          i0 = loop_no*n_grid_points_per_loop + j*n_grid_points_x+i
          i1 = loop_no*n_grid_points_per_loop + j*n_grid_points_x+i+1
          i2 = loop_no*n_grid_points_per_loop + j*n_grid_points_x+i+2
          
          i3 = loop_no*n_grid_points_per_loop + (j+1)*n_grid_points_x+i
          i4 = loop_no*n_grid_points_per_loop + (j+1)*n_grid_points_x+i+1
          i5 = loop_no*n_grid_points_per_loop + (j+1)*n_grid_points_x+i+2
          
          i6 = loop_no*n_grid_points_per_loop + (j+2)*n_grid_points_x+i
          i7 = loop_no*n_grid_points_per_loop + (j+2)*n_grid_points_x+i+1
          i8 = loop_no*n_grid_points_per_loop + (j+2)*n_grid_points_x+i+2
          
          i9 = (loop_no+1)*n_grid_points_per_loop + j*n_grid_points_x+i
          i10 = (loop_no+1)*n_grid_points_per_loop + j*n_grid_points_x+i+1
          i11 = (loop_no+1)*n_grid_points_per_loop + j*n_grid_points_x+i+2
          
          i12 = (loop_no+1)*n_grid_points_per_loop + (j+1)*n_grid_points_x+i
          i13 = (loop_no+1)*n_grid_points_per_loop + (j+1)*n_grid_points_x+i+1
          i14 = (loop_no+1)*n_grid_points_per_loop + (j+1)*n_grid_points_x+i+2
          
          i15 = (loop_no+1)*n_grid_points_per_loop + (j+2)*n_grid_points_x+i
          i16 = (loop_no+1)*n_grid_points_per_loop + (j+2)*n_grid_points_x+i+1
          i17 = (loop_no+1)*n_grid_points_per_loop + (j+2)*n_grid_points_x+i+2
         
          i18 = (loop_no+2)*n_grid_points_per_loop + j*n_grid_points_x+i
          i19 = (loop_no+2)*n_grid_points_per_loop + j*n_grid_points_x+i+1
          i20 = (loop_no+2)*n_grid_points_per_loop + j*n_grid_points_x+i+2
          
          i21 = (loop_no+2)*n_grid_points_per_loop + (j+1)*n_grid_points_x+i
          i22 = (loop_no+2)*n_grid_points_per_loop + (j+1)*n_grid_points_x+i+1
          i23 = (loop_no+2)*n_grid_points_per_loop + (j+1)*n_grid_points_x+i+2
          
          i24 = (loop_no+2)*n_grid_points_per_loop + (j+2)*n_grid_points_x+i
          i25 = (loop_no+2)*n_grid_points_per_loop + (j+2)*n_grid_points_x+i+1
          i26 = (loop_no+2)*n_grid_points_per_loop + (j+2)*n_grid_points_x+i+2
         
          quadratic_elements.append([i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26])
    
  if debugging_stl_output:
    # create stl mesh of 3D mesh
    for linear_element in linear_elements:
      
      p = []
      for i in range(8):
        p.append([])
        p[i] = nodes[linear_element[i]]
      
      out_3d_mesh_triangles += [
        [p[0],p[3],p[1]],[p[0],p[2],p[3]],  # bottom
        [p[4],p[5],p[7]],[p[4],p[7],p[6]],  # top
        [p[0],p[1],p[5]],[p[0],p[5],p[4]],  # front
        [p[2],p[7],p[3]],[p[2],p[6],p[7]],  # back
        [p[2],p[0],p[4]],[p[2],p[4],p[6]],  # left
        [p[1],p[3],p[7]],[p[1],p[7],p[5]]  # right
      ]
    
  # create seed points as the 2D linear elements' center points in the middle loop
  seed_points = []
  middle_loop_no = int(n_loops/2)
  middle_loop_grid_points = loop_grid_points[middle_loop_no]

  # linear elements
  for j in range(n_grid_points_y-1):
    for i in range(n_grid_points_x-1):
      
      i0 = j*n_grid_points_x+i
      i1 = j*n_grid_points_x+i+1
      i2 = (j+1)*n_grid_points_x+i
      i3 = (j+1)*n_grid_points_x+i+1
      
      center_point = (middle_loop_grid_points[i0] + middle_loop_grid_points[i1] + middle_loop_grid_points[i2] + middle_loop_grid_points[i3])/4.
      seed_points.append(center_point.tolist())
    
  # compute indices of bottom and top nodes. These are needed for setting boundary conditions
  bottom_node_indices = list(range(n_grid_points_per_loop))
  top_node_indices = list(range(n_grid_points_per_loop*(n_loops-1),n_grid_points_per_loop*n_loops))
    
  # set number of elements per coordinate direction, which is needed for structured meshes
  n_linear_elements_per_coordinate_direction = [n_grid_points_x-1, n_grid_points_y-1, n_loops-1]
  n_quadratic_elements_per_coordinate_direction = [(n_grid_points_x-1)/2, (n_grid_points_y-1)/2, (n_loops-1)/2]
    
  # output nodes and elements and seed points
  data = {
    "node_positions": node_positions, 
    "linear_elements": linear_elements, 
    "quadratic_elements": quadratic_elements, 
    "seed_points": seed_points,
    "bottom_nodes": bottom_node_indices,
    "top_nodes": top_node_indices,
    "n_linear_elements_per_coordinate_direction": n_linear_elements_per_coordinate_direction,
    "n_quadratic_elements_per_coordinate_direction": n_quadratic_elements_per_coordinate_direction,
  }
  return data
  
def create_3d_mesh_simple(n_points_x, loops):
  """
  This performs the whole process of generating the 3D mesh, starting from the rings
  :param n_points_x: number of points per dimension for the initial triangulation, assuming a cube reference area (4*n_points_x on circumference for unit circle)
  """
    
  # constant parameters
  triangulation_type = 2  # 0 = scipy, 1 = triangle, 2 = center pie (2 is best), 3 = minimized distance
  parametric_space_shape = 3   # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
  max_area_factor = 2.    # only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
  show_plot = False
  debug = False  
  debugging_stl_output = False
  improve_mesh = False

  n_grid_points_x = n_points_x+1   # grid width of generated 2d mesh
  n_grid_points_y = n_points_x+1

  if parametric_space_shape == 0:  # for unit circle 
    n_grid_points_y = 20
        
  n_loops = len(loops)
  print("{} loops".format(n_loops))

  # sample loop with 4*n_points_x equidistant points
  n_points = 4*n_points_x
  border_point_loops,lengths = rings_to_border_points(loops, n_points)

  loop_grid_points = []  # list of grid point, for every slice, only contains loops that are not empty  
    
  # loop over all loops of border points
  for loop_no,(border_points,length) in enumerate(zip(border_point_loops,lengths)):
    
    print("")
    print("Ring {}/{} with {} border points, length: {}".format(loop_no, n_loops, len(border_points), length))
    
    # create 2D mesh with border_points
    show_plot = False
    grid_points_world_space,duration_1d = stl_create_mesh.create_planar_mesh(border_points, loop_no, n_points, \
      n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, improve_mesh, show_plot, debugging_stl_output, [])
        
    # store grid points in world space of current loop
    loop_grid_points.append(grid_points_world_space)

  data = stl_create_mesh.create_3d_mesh(loop_grid_points, n_grid_points_x, n_grid_points_y, debugging_stl_output, [])
  return data

def create_3d_mesh_from_border_points_faces(border_points_faces, improve_mesh, level_no):
  """
  Create the 3D mesh from border points which are organised as faces.
  :param border_points_faces: [border_points_0minus, border_points_0plus, border_points_1minus, border_points_1plus]
  :param improve_mesh: if the 2D meshes should be smoothed, this takes a lot of time but improves the result
  """
  
  print("create_3d_mesh_from_border_points_faces")

  if False:
    filename = "dump_border_points_faces_{}.py".format(os.getpid())
    print("dump to filename:{}".format(filename))
    with open(filename, 'wb') as f:
      pickle.dump(border_points_faces, f)
          
  # constant parameters
  triangulation_type = 1  # 0 = scipy, 1 = triangle, 2 = center pie (2 is best), 3 = minimized distance
  parametric_space_shape = 3   # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
  max_area_factor = 2.    # only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
  show_plot = False
  debugging_stl_output = False
  #improve_mesh = True    # post-smooth mesh

  border_points_0minus = border_points_faces[0]   # the first / last point of each list for the face overlaps with an identical point on another face's list
  border_points_0plus = border_points_faces[1]
  border_points_1minus = border_points_faces[2]
  border_points_1plus = border_points_faces[3]
  n_loops = len(border_points_0minus)
  #print("border_points_0minus: {}".format(border_points_0minus))
  
  #   ^ --(1+)-> ^
  # ^ 0-         0+
  # | | --(1-)-> |
  # +-->

  border_point_loops = []
  for loop_no in range(n_loops):
    border_point_loop = []
    
    # do not consider last point of border_points which is the same as the first of the next face
    for point in border_points_1minus[loop_no][0:-1]:
      border_point_loop.append(np.array(point))

    for point in border_points_0plus[loop_no][0:-1]:
      border_point_loop.append(np.array(point))

    for i,point in enumerate(reversed(border_points_1plus[loop_no])):
      if i == len(border_points_1plus[loop_no])-1:
        break
      border_point_loop.append(np.array(point))

    for i,point in enumerate(reversed(border_points_0minus[loop_no])):
      if i == len(border_points_0minus[loop_no])-1:
        break
      border_point_loop.append(np.array(point))

    border_point_loops.append(border_point_loop)

  #print("border_point_loops: {}".format(border_point_loops))

  n_points = len(border_point_loops[0])    
  n_points_x = (int)(n_points/4)  
  n_grid_points_x = n_points_x+1   # grid width of generated 2d mesh
  n_grid_points_y = n_points_x+1
    
  def handle_loop(loop_no, border_points):
    n_points = len(border_points)    
    n_points_x = (int)(n_points/4)  
    n_grid_points_x = n_points_x+1   # grid width of generated 2d mesh
    n_grid_points_y = n_points_x+1

    print("")
    print("Level {}, ring {}/{} with {} border points".format(level_no, loop_no, n_loops, n_points))
    #print("border points: ", border_points)
    
    # create 2D mesh with border_points
    if debugging_stl_output:
      out_triangulation_world_space = []
      markers_border_points_world_space = []
      out_triangulation_parametric_space = []
      grid_triangles_world_space = []
      grid_triangles_parametric_space = []
      markers_grid_points_parametric_space = []
      markers_grid_points_world_space = []
      debugging_output_lists = [out_triangulation_world_space, markers_border_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,markers_grid_points_parametric_space, markers_grid_points_world_space]
    else:
      debugging_output_lists = []
    
    if loop_no == 0 and False:
      output = (border_points, loop_no, n_points, n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, show_plot)
      filename = "dump_mesh_2d_{}.py".format(os.getpid())
      print("dump to filename:{}".format(filename))
      with open(filename, 'wb') as f:
        pickle.dump(output, f)
    
    grid_points_world_space,duration_1d = create_planar_mesh(border_points, loop_no, n_points, \
      n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, improve_mesh, show_plot, debugging_stl_output, debugging_output_lists)
      
    if debugging_stl_output and False:
      stl_debug_output.output_triangles("2dmesh_loop_{}_w_triangulation".format(loop_no), out_triangulation_world_space)
      stl_debug_output.output_triangles("2dmesh_loop_{}_p_triangulation".format(loop_no), out_triangulation_parametric_space)
      stl_debug_output.output_triangles("2dmesh_loop_{}_p_grid".format(loop_no), grid_triangles_parametric_space + markers_grid_points_parametric_space)
      stl_debug_output.output_triangles("2dmesh_loop_{}_w_border".format(loop_no), markers_border_points_world_space)
      stl_debug_output.output_triangles("2dmesh_loop_{}_w_grid".format(loop_no), grid_triangles_world_space + markers_grid_points_world_space)
        
    # store grid points in world space of current loop
    #loop_grid_points.append(grid_points_world_space)
    return grid_points_world_space
  
  loop_grid_points = [[] for i in range(len(border_point_loops))]  # list of grid point, for every slice, only contains loops that are not empty  
  
  print("create 2D meshes on {} loops, then create a 3D mesh".format(len(border_point_loops)))
  # serial implementation
  if True:
    for loop_no,border_points in enumerate(border_point_loops):
      loop_grid_points[loop_no] = handle_loop(loop_no, border_points)
      
      if len(loop_grid_points[loop_no]) == 0:
        return None
  
  # concurrent execution, currently not working, but not clear why
  if False:
    print("start concurrent call")
    # concurrently call handle_loop(loop_no, border_points) with the points in border_point_loops, this is the same as 
    # for loop_no,border_points in enumerate(border_point_loops):
    #   loop_grid_points[loop_no] = handle_loop(loop_no, border_points)
    #
    try:
      import concurrent.futures
      with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
        
        print("create futures")
        futures = {executor.submit(handle_loop, loop_no, border_points): loop_no for loop_no,border_points in enumerate(border_point_loops)}
        print("futures created")
        #concurrent.futures.wait(futures)
        #print("futures finished")
        for future in concurrent.futures.as_completed(futures):
          loop_no = futures[future]
          print("{} has completed".format(loop_no))
          
          try:
            result = future.result()
          except Exception as exc:
            print("Ring {} generated an exception: {}".format(loop_no, exc))
          else:
            print("Ring {}/{} finished".format(loop_no,n_loops))
            loop_grid_points[loop_no] = result
    except Exception as e:
      print(e)
    else:
      print("success")
    print("after concurrent futures")
  
  triangles_3dmesh = []
  data = create_3d_mesh(loop_grid_points, n_grid_points_x, n_grid_points_y, debugging_stl_output, triangles_3dmesh)
  
  if debugging_stl_output:
    stl_debug_output.output_triangles("3dmesh", triangles_3dmesh)
  
  return data
