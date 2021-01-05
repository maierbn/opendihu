#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# This file contains functionality to creates a hexaeder mesh out of the rings given by create_rings.py. 
# The terms "ring" and "loop" are used synonymously in this file.
# Use create_mesh.py to actual perform the whole pipeline, or look at create_3d_mesh_simple or create_3d_mesh_from_boundary_points_faces.
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
import create_planar_mesh_helper
# helper functions for create_planar_mesh
import create_slice_tringulation
import solve_laplace_problem
import map_quadrangulation_to_world_space
import fix_and_smooth_mesh

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
    
def sample_boundary_points(loop, length, n_points, target_x, last_loop_start_point):
  """
  Given a loops of points (which are spaced like the initial STL triangles were), sample the loops at n_points equidistant points.
  The starting point is determined as follows: If last_loop_start_point is None, it determines the start point from target_x.
  The start point lies on the intersection of the plane x=target_x with the loop. if there are multiple intersections, the first that is in y- direction from the cog is taken.
  If last_loop_start_point is given, the start point is the point which has the lowest distance to last_loop_start_point., i.e. the connection line is perpendicular to the loop.
  This is helpful to ensure that two neighbouring loops on different z-levels yield similar boundary point loops.
  
  :param loop:   List of points that form the loop to be sampled.
  :param length:  The length of the loop, should be measured previously
  :param n_points: Number of points that will be created.
  :param target_x: a plane x=target_x that will determine the start point in the loop. The start point lies on that plane and is the one left (in x- direction) of the center.
  :param last_loop_start_point: If None, target_x will be used to determine the start_point, if given the start_point is chosen in a way such that its distance to last_loop_start_point is minimal.
  :return: start_point, boundary_points
  """
  
  debug = False
  boundary_points = []            # this list will hold the result
  h = float(length) / n_points  # parameter increment between boundary points
  
  # determine start point
  center_point = sum(loop)/len(loop)
  if target_x is None:
    target_x = 0.5*(loop[0][0]+loop[-1][0])   # plane x=target_x on which the start point will lie, it has to be ensured, that all rings touch this plane
  
  if debug:
    print("sample_boundary_points, loop with {} points, length: {}, sample {} points, target_x: {}, h: {}, last_loop_start_point:{}\n".format(len(loop), length, n_points, target_x, h, last_loop_start_point))
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
  
  # iterate over the loop of points until the start point is reached, then iterate one more time, collecting the sampled boundary points in equidistant intervals
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
        # now the search for the start point has ended and the normal iterations where boundary points will be sampled continues
        start_point_reached = True
        
        # if the start point is not yet set
        if start_point is None:
          edge = -previous_loop_point + loop_point 
          edge_length = np.linalg.norm(edge)
          alpha = (target_x - previous_loop_point[0]) / edge_length
          start_point = previous_loop_point + alpha*edge
        
        # add starting point of loop as first boundary points
        boundary_points.append(start_point)
        if debug:
          print("   add start point: ",boundary_points)
        
        previous_loop_point = start_point
        
        t_previous_loop_point = 0
        t_next_boundary_point = h
          
    # if we are at the second part of the loop iteration where the actual boundary points are collected
    if start_point_reached:
      # compute current edge
      edge = -previous_loop_point + loop_point 
      edge_length = np.linalg.norm(edge) 
      
      n_on_edge = 0
      
      # collect all boundary points on current edge
      while t_next_boundary_point <= t_previous_loop_point + edge_length:
        
        boundary_point = previous_loop_point + edge * (t_next_boundary_point - t_previous_loop_point) / edge_length
        boundary_points.append(boundary_point)
        
        if debug:
          print("   > boundary point {}, t previous: {}, edge length: {}, next: {},  t_previous_loop_point + edge_length: {}".format(boundary_point, t_previous_loop_point, edge_length, t_next_boundary_point,  t_previous_loop_point + edge_length))
        
        t_next_boundary_point += h
        if debug:
          n_on_edge += 1
          print("   new t_next_boundary_point: ",t_next_boundary_point)
      
      if debug:
        print("   edge length {}, edge ({}, {}), n_on_edge: {}, t_previous_loop_point: {} t_next_boundary_point: {}".format(edge_length,previous_loop_point,loop_point,n_on_edge,t_previous_loop_point,t_next_boundary_point))
      
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
      return sample_boundary_points(loop, length, n_points, target_x, last_loop_start_point)
  
  # end of while loop
  
  if debug:
    print("visited ring length: {}, actual length: {}".format(t_previous_loop_point, length))
  
  if debug:
    print("distances between boundary points, note this are not necessary the same as h:")
    p0 = boundary_points[0]
    for p in boundary_points[1:]:
      print("  {}".format(np.linalg.norm(p0-p)))
      p0 = p
  
  # If there were too many points collected, crop to needed points.
  # This happens because on the edge of the starting point, p0-ps-p1 (with edge p0-p1 and ps=starting point), the part ps-p1 is visited twice.
  if len(boundary_points) > n_points:
    if debug:
      print("too many points: {}, n_points: {}, {}\n\n".format(len(boundary_points), n_points, boundary_points))
    boundary_points = boundary_points[:n_points]
  elif len(boundary_points) < n_points:
    print("Error, not enough points: {}, n_points: {}, {}\n\n".format(len(boundary_points), n_points, boundary_points))
  
  if debug:
    print("result: ")
    print("start point: ",start_point)
    print("boundary points: ",len(boundary_points))
    print(boundary_points)
    print("h: ",h)
    print("n boundary points: {}".format(len(boundary_points)))
  
  return start_point,boundary_points
    
def rings_to_boundary_points(loops, n_points):
  """
  Standardize every ring to be in counter-clockwise direction and starting with the point with lowest x coordinate,
  then sample boundary points
  :param loops: list of rings
  :param n_points: number of boundary points that should be sampled per ring
  :returns: boundary_point_loops,lengths  where lengths is a list of lengths of each loop and boundary_point_loops is the list of loops with boundary points
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
  
  # sample boundary points for each ring
  boundary_point_loops = []
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
    start_point,boundary_points = sample_boundary_points(loop, length, n_points, target_x, last_loop_start_point)
    last_loop_start_point = start_point
    boundary_point_loops.append(boundary_points)
    
  return boundary_point_loops,lengths

def boundary_point_loops_to_list(boundary_point_loops):
  """
  transform the points from numpy array to list, such that they can be extracted from the opendihu C++ code
  """
  result = []
  for boundary_point_loop in boundary_point_loops:
    result_boundary_points = []
    for boundary_point in boundary_point_loop:
      result_boundary_points.append(boundary_point.tolist())
    result.append(result_boundary_points)
  return result

def create_planar_mesh(boundary_points, loop_no, n_points, \
    n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, improve_mesh, show_plot, debugging_stl_output, stl_triangle_lists):
  """
  Create a 2d mesh using the harmonic map algorithm. Input is a list of boundary points which will be the boundary nodes of the mesh.
  :param boundary_points: list of points on the boundary that will be the nodes, (each point is again a list with 3 entries)
  :param loop_no: the number of the current loop, this is only needed for debugging output filename
  :param n_points: number of points on the boundary of the resulting mesh
  :param n_grid_points_x: number of grid points in x direction of the final quadrangulation
  :param n_grid_points_y: number of grid points in y direction of the final quadrangulation
  :param triangulation_type: 0 = scipy, 1 = triangle, 2 = center pie (2 is best), 3 = minimized distance
  :param parametric_space_shape: 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
  :param max_area_factor: only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
  :param improve_mesh: if the generated mesh should be smoothed afterwards
  :param show_plot: if the matplotlib plots of geometry and harmonic map should be generated
  :param debugging_stl_output: if list should be filled with STL triangles that can be output to a STL mesh for debugging
  :param stl_triangle_lists: the debugging lists: [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,markers_grid_points_parametric_space, markers_grid_points_world_space]
  """

  # note: the first and last point in the loop is the same, because it represents the edges between the points

  global distances_between_world_mesh_nodes_std   # list of distances between neighboring grid points
  global relative_distances_between_world_mesh_nodes_std   # list of relative distances between neighbouring grid points, relative per slice
  
  t_start = timeit.default_timer()

  debug = False   # enable debugging output
  
  if debugging_stl_output:
    [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,\
      markers_grid_points_parametric_space, markers_grid_points_world_space] = stl_triangle_lists
    
  if len(boundary_points) != n_points:
    print("\nError in create_planar_mesh, n_points: {}, but {} boundary points given.\n boundary points: {}\n".format(n_points, len(boundary_points), boundary_points))
  
  modify_phi = False
  if parametric_space_shape == 3:
    modify_phi = True
  
  # triangulate surface in world space
  points = np.reshape(boundary_points,(n_points,3))
  n_points_per_face = (int)(n_points/4)
  
  if debug:  
    print("")
    print("points:")
    print(points)

  if debugging_stl_output:
    # for debugging create markers at boundary points
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
      markers_boundary_points_world_space += [
        [point+diag0,point+diag3,point+diag1],[point+diag0,point+diag2,point+diag3],  # bottom
        [point+diag4,point+diag5,point+diag7],[point+diag4,point+diag7,point+diag6],  # top
        [point+diag0,point+diag1,point+diag5],[point+diag0,point+diag5,point+diag4],  # front
        [point+diag2,point+diag7,point+diag3],[point+diag2,point+diag6,point+diag7],  # back
        [point+diag2,point+diag0,point+diag4],[point+diag2,point+diag4,point+diag6],  # left
        [point+diag1,point+diag3,point+diag7],[point+diag1,point+diag7,point+diag5]  # right
      ]

  # In the following, helper functions are called that implement the high-level steps of the algorithm. These helper functions are defined in their own script files.
  # The high number of return variables is because the output_debugging_files functions needs them to output them all

  # create a triangulation on the 2D slice
  point_indices_list, triangle_list, n_points, points, n_original_points, original_points, original_point_phi_value, get_modified_phi, n_regular_grid_boundary_points, extent_x, extent_y, n_additional_points_on_ring, determine_additional_points_on_ring \
   = create_slice_tringulation.create_slice_tringulation(triangulation_type, modify_phi, n_points, points, max_area_factor)

  # solve Laplace equation on deformed mesh
  u,v = solve_laplace_problem.solve_laplace_problem(n_points, points, n_original_points, original_points, original_point_phi_value, triangulation_type, parametric_space_shape, point_indices_list, triangle_list, debugging_stl_output, stl_triangle_lists)
  
  # now the mapping x -> u,v is computed (from parameter space to world space)
  grid_points_world_space, grid_points_parametric_space, grid_points_parametric_space_modified \
   = map_quadrangulation_to_world_space.map_quadrangulation_to_world_space(u, v, n_grid_points_x, n_grid_points_y, parametric_space_shape, modify_phi, original_point_phi_value, get_modified_phi, n_regular_grid_boundary_points, point_indices_list, triangle_list, extent_x, points, debugging_stl_output, stl_triangle_lists)
  
  # improve the world mesh by Laplacian smoothing and fixing of self-intersections
  grid_points_world_space_improved = grid_points_world_space
  
  # if this is enabled
  if improve_mesh:   
    # optimize points, fix self-intersections and Laplacian smoothing
    grid_points_world_space_improved \
     = fix_and_smooth_mesh.fix_and_smooth_mesh(grid_points_world_space, n_grid_points_x, n_grid_points_y, point_indices_list, triangle_list, extent_x, extent_y, loop_no, debugging_stl_output, stl_triangle_lists)
  
  # create various STL files that help understand and debug the algorithm
  if debugging_stl_output:
    create_planar_mesh_helper.output_debugging_files(grid_points_parametric_space, grid_points_world_space, grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, parametric_space_shape, 
                                                     u, v, points, n_points_per_face, n_additional_points_on_ring, loop_no, show_plot, extent_x, determine_additional_points_on_ring, grid_points_parametric_space_modified, point_indices_list, triangle_list, stl_triangle_lists)
  
  t_stop = timeit.default_timer()
  duration = t_stop - t_start
  
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
  boundary_point_loops,lengths = rings_to_boundary_points(loops, n_points)

  loop_grid_points = []  # list of grid point, for every slice, only contains loops that are not empty  
    
  # loop over all loops of boundary points
  for loop_no,(boundary_points,length) in enumerate(zip(boundary_point_loops,lengths)):
    
    print("")
    print("Ring {}/{} with {} boundary points, length: {}".format(loop_no, n_loops, len(boundary_points), length))
    
    # create 2D mesh with boundary_points
    show_plot = False
    grid_points_world_space,duration_1d = stl_create_mesh.create_planar_mesh(boundary_points, loop_no, n_points, \
      n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, improve_mesh, show_plot, debugging_stl_output, [])
        
    # store grid points in world space of current loop
    loop_grid_points.append(grid_points_world_space)

  data = stl_create_mesh.create_3d_mesh(loop_grid_points, n_grid_points_x, n_grid_points_y, debugging_stl_output, [])
  return data

def create_3d_mesh_from_boundary_points_faces(boundary_points_faces, improve_mesh, max_area_factor, level_no):
  """
  Create the 3D mesh from boundary points which are organised as faces.
  :param boundary_points_faces: [boundary_points_0minus, boundary_points_0plus, boundary_points_1minus, boundary_points_1plus]
  :param improve_mesh: if the 2D meshes should be smoothed, this takes a lot of time but improves the result
  :param max_area_factor: only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
  :param level_no: only for debugging output
  """
  
  if False:
    filename = "dump_boundary_points_faces_{}.py".format(os.getpid())
    print("dump to filename:{}".format(filename))
    with open(filename, 'wb') as f:
      pickle.dump(boundary_points_faces, f)
          
  # constant parameters
  triangulation_type = 2  # 0 = scipy, 1 = triangle, 2 = center pie (2 is best), 3 = minimized distance
  parametric_space_shape = 3   # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
  #max_area_factor = 100.    # only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
  show_plot = False               # set this to true such that the plots will be opened instead of written to a file
  debugging_stl_output = False       # set this value to true to enable the matplotlib output during the parallel_fiber_estimation runs
  #improve_mesh = True    # post-smooth mesh

  boundary_points_0minus = boundary_points_faces[0]   # the first / last point of each list for the face overlaps with an identical point on another face's list
  boundary_points_0plus = boundary_points_faces[1]
  boundary_points_1minus = boundary_points_faces[2]
  boundary_points_1plus = boundary_points_faces[3]
  n_loops = len(boundary_points_0minus)
  #print("boundary_points_0minus: {}".format(boundary_points_0minus))
  
  #   ^ --(1+)-> ^
  # ^ 0-         0+
  # | | --(1-)-> |
  # +-->

  boundary_point_loops = []
  for loop_no in range(n_loops):
    boundary_point_loop = []
    
    # do not consider last point of boundary_points which is the same as the first of the next face
    for point in boundary_points_1minus[loop_no][0:-1]:
      boundary_point_loop.append(np.array(point))

    for point in boundary_points_0plus[loop_no][0:-1]:
      boundary_point_loop.append(np.array(point))

    for i,point in enumerate(reversed(boundary_points_1plus[loop_no])):
      if i == len(boundary_points_1plus[loop_no])-1:
        break
      boundary_point_loop.append(np.array(point))

    for i,point in enumerate(reversed(boundary_points_0minus[loop_no])):
      if i == len(boundary_points_0minus[loop_no])-1:
        break
      boundary_point_loop.append(np.array(point))

    boundary_point_loops.append(boundary_point_loop)

  #print("boundary_point_loops: {}".format(boundary_point_loops))

  n_points = len(boundary_point_loops[0])    
  n_points_x = (int)(n_points/4)  
  n_grid_points_x = n_points_x+1   # grid width of generated 2d mesh
  n_grid_points_y = n_points_x+1
    
  def handle_loop(loop_no, boundary_points):
    n_points = len(boundary_points)    
    n_points_x = (int)(n_points/4)  
    n_grid_points_x = n_points_x+1   # grid width of generated 2d mesh
    n_grid_points_y = n_points_x+1

    #print("")
    print("Level {}, ring {}/{} with {} boundary points".format(level_no, loop_no, n_loops, n_points))
    #print("boundary points: ", boundary_points)
    
    # create 2D mesh with boundary_points
    if debugging_stl_output:
      out_triangulation_world_space = []
      markers_boundary_points_world_space = []
      out_triangulation_parametric_space = []
      grid_triangles_world_space = []
      grid_triangles_parametric_space = []
      markers_grid_points_parametric_space = []
      markers_grid_points_world_space = []
      debugging_output_lists = [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,markers_grid_points_parametric_space, markers_grid_points_world_space]
    else:
      debugging_output_lists = []
    
    if loop_no == 0 and False:
      output = (boundary_points, loop_no, n_points, n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, show_plot)
      filename = "dump_mesh_2d_{}.py".format(os.getpid())
      print("dump to filename:{}".format(filename))
      with open(filename, 'wb') as f:
        pickle.dump(output, f)
    
    grid_points_world_space,duration_1d = create_planar_mesh(boundary_points, loop_no, n_points, \
      n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, improve_mesh, show_plot, debugging_stl_output, debugging_output_lists)
      
    if debugging_stl_output and False:
      stl_debug_output.output_triangles("2dmesh_loop_{}_w_triangulation".format(loop_no), out_triangulation_world_space)
      stl_debug_output.output_triangles("2dmesh_loop_{}_p_triangulation".format(loop_no), out_triangulation_parametric_space)
      stl_debug_output.output_triangles("2dmesh_loop_{}_p_grid".format(loop_no), grid_triangles_parametric_space + markers_grid_points_parametric_space)
      stl_debug_output.output_triangles("2dmesh_loop_{}_w_boundary".format(loop_no), markers_boundary_points_world_space)
      stl_debug_output.output_triangles("2dmesh_loop_{}_w_grid".format(loop_no), grid_triangles_world_space + markers_grid_points_world_space)
        
    # store grid points in world space of current loop
    #loop_grid_points.append(grid_points_world_space)
    return grid_points_world_space
  
  loop_grid_points = [[] for i in range(len(boundary_point_loops))]  # list of grid point, for every slice, only contains loops that are not empty  
  
  print("create 2D meshes on {} loops, then create a 3D mesh".format(len(boundary_point_loops)))
  # serial implementation
  if True:
    for loop_no,boundary_points in enumerate(boundary_point_loops):
      loop_grid_points[loop_no] = handle_loop(loop_no, boundary_points)
      
      if len(loop_grid_points[loop_no]) == 0:
        return None
  
  # concurrent execution, currently not working, but not clear why
  if False:
    print("start concurrent call")
    # concurrently call handle_loop(loop_no, boundary_points) with the points in boundary_point_loops, this is the same as 
    # for loop_no,boundary_points in enumerate(boundary_point_loops):
    #   loop_grid_points[loop_no] = handle_loop(loop_no, boundary_points)
    #
    try:
      import concurrent.futures
      with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
        
        print("create futures")
        futures = {executor.submit(handle_loop, loop_no, boundary_points): loop_no for loop_no,boundary_points in enumerate(boundary_point_loops)}
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
