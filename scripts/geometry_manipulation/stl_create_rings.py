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
import stl_create_mesh   # for standardize_loop and rings_to_boundary_points
import stl_debug_output
import spline_surface

import stl
from stl import mesh
from svg.path import parse_path
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier

def get_intersecting_line_segment(triangle, z_value):
  """ Intersect a triangle with a z plane. 
      For a triangle T given by points p1,p2,p3 ∈ ℝ^3 and a value z, determine the the set of points P = T ∩ {p | p_z=z},
      i.e., the line segment [pa, pb] that lies in the triangle T and on the z plane z=z_value """

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
  
  # if one corner point lies on the z plane
  if abs(p1z-z_value) < 1e-12:
    return [p1,p1]
  if abs(p2z-z_value) < 1e-12:
    return [p2,p2]
  if abs(p3z-z_value) < 1e-12:
    return [p3,p3]

  # if p2z == p1z swap p1 and p3 
  if p2z == p1z:
    t = p1
    p1 = p3
    p3 = t
  
    p1z = p1[2]
    p3z = p3[2]
  
  debug = False
  
  # for debugging, determine if the triangle has to intersect the plane, i.e. there has to be an intersection point
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
  
  # check which of the boundarys of the triangle in parameter space are intersected by the line segment intersects
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

def create_loop(z_value, stl_mesh, loop):
  """
  Sample the mesh at a given z level, result is a loop, the edges are in no particular order.
  The output loop is a list of edges.
  :param z_value: z level of all extracted points
  :param stl_mesh: the mesh as stl.mesh object
  :param loop: this is the output, this list will contain edges [p0,p1]
  """
  debug = False

  if debug:
    print(" z_value: {}, n points: {}: {}".format(z_value, len(stl_mesh.points), stl_mesh.points))

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
      for [p1, p2] in loop:
        if (np.allclose(p1, edge[0]) and np.allclose(p2, edge[1])) \
          or (np.allclose(p2, edge[0]) and np.allclose(p1, edge[1])):
          edge_is_already_in_loop = True
          break
      
      
      # append edge to loop
      if not edge_is_already_in_loop:
        if debug:
          print(" add edge ",edge," to loop")
        #print(", prev: ", loop, "->", 
        loop.append(edge)
        #print(loop
      else:
        if debug:
          print("edge_is_already_in_loop")
      
def order_loop(loop, first_point):
  """
  Create a new loop with consecutive edges, starting from first_point, which has to be a point contained in the loop
  :param loop: list of edges in no particular order
  :param first_point: the point which will be the first point of the resulting loop, must be a point of the input loop
  :return: a list of points
  """
  
  debug = False
  
  # fill new loop starting with first_point, ordered by consecutive edges
  new_loop = [first_point]
  
  if debug:
    print("")
    #print("order_loop(loop={},first_point={})".format(loop,first_point))
    print("first point: ", first_point)
  
  points_with_single_adjacent_edge = []
  
  previous_end_point = first_point
  current_end_point = first_point
  
  while len(new_loop) < len(loop)+2:
      
    next_point_found = False
      
    # iterate over points in old loop
    for edge in loop:
      
      if np.allclose(current_end_point, edge[0]) and not np.allclose(previous_end_point, edge[1]):
        # add edge p0->p1 to new_loop
        new_loop.append(edge[1])
        previous_end_point = current_end_point
        current_end_point = edge[1]
        next_point_found = True
        
        if debug:
          print("add point ",edge[1],", new length of loop: ",len(new_loop),", expected final length: ",len(loop)+1)
          
      elif np.allclose(current_end_point, edge[1]) and not np.allclose(previous_end_point, edge[0]):
        # add edge p1->p0 to new_loop
        new_loop.append(edge[0])
        previous_end_point = current_end_point
        current_end_point = edge[0]
        next_point_found = True
        
        if debug:
          print("add point ",edge[0],", new length of loop: ",len(new_loop),", expected final length: ",len(loop)+1)
        
    if not next_point_found:
      if debug: 
        print("no point found that continues loop")
      
      # detect points that have only one adjacent edge, if not already done
      if len(points_with_single_adjacent_edge) == 0:
        
        # iterate over points in old loop
        for edge in loop:
          for p in edge:
            n_adjacent_edges = 0
            # iterate over all edges
            for edge2 in loop:
              if np.allclose(p, edge2[0]) or np.allclose(p, edge2[1]):
                n_adjacent_edges += 1
              if n_adjacent_edges >= 2:
                break
            if n_adjacent_edges == 1:
              points_with_single_adjacent_edge.append(p)
      
      # now, points_with_single_adjacent_edge contains all points that have only a single adjacent edge
      
      # look for closest such point and use that as next one
      minimum_distance = None
      
      closest_point = None
      # iterate over points in points_with_single_adjacent_edge
      for point in points_with_single_adjacent_edge:
        # do not allow last edge to be selected again
        if not np.allclose(current_end_point, point) and not np.allclose(previous_end_point, point) :
          distance = np.linalg.norm(point - current_end_point)
          if minimum_distance is None or distance < minimum_distance:
            minimum_distance = distance
            closest_point = point
          
      new_loop.append(closest_point)
      previous_end_point = current_end_point
      current_end_point = closest_point
      if debug: 
        print("use closest point {} with a distance of {} to the current point".format(closest_point, minimum_distance))
       
    # if the end point is again the start point, finish the loop. Note if by now still len(new_loop) < len(loop)+1 holds, there might be a different (distinct) loop for this z position which is discarded.
    if np.allclose(current_end_point, first_point):
      if debug:
        print("start point reached")
      break
       
      #print("Error: Loop for z={} could not be closed. Maybe there are triangles missing?".format(loop[0][0][2]))
      #break
  return new_loop
      
def create_rings(input_filename, bottom_clip, top_clip, n_loops, write_output_mesh):
  """
  Create n_loops rings/loops (slices) on a closed surface, in equidistant z-values between bottom_clip and top_clip (including those)
  :param input_filename: file name of an stl file that contains the closed surface mesh of the muscle, aligned with the z-axis
  :param bottom_clip: the bottom z-value where to clip the muscle
  :param top_clip: the top z-value where to clip the muscle
  :param n_loops: number of loops/rings to extract between bottom_clip and top_clip
  :param write_output_mesh: if the result should be written to an output stl file for debugging
  :return: a list of lists of points (one list for each level)
  """

  stl_mesh = mesh.Mesh.from_file(input_filename)

  out_triangles = []

  n_is_inside_1 = 0
  n_is_inside_2 = 0

  n_triangles = len(stl_mesh.points)
  
  # determine bounding box of mesh
  bounding_box_zmin = np.inf
  bounding_box_zmax = -np.inf
  # loop over triangles in mesh
  for p in stl_mesh.points:
    # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

    p1 = np.array(p[0:3])
    p2 = np.array(p[3:6])
    p3 = np.array(p[6:9])
    
    bounding_box_zmin = min(bounding_box_zmin,p1[2],p2[2],p3[2])
    bounding_box_zmax = max(bounding_box_zmax,p1[2],p2[2],p3[2])
  
  if bottom_clip < bounding_box_zmin:
    print("Adjusting bottom_clip from {} to {} to match the mesh.".format(bottom_clip,bounding_box_zmin))
    bottom_clip = bounding_box_zmin
  if top_clip > bounding_box_zmax:
    print("Adjusting top_clip from {} to {} to match the mesh.".format(top_clip,bounding_box_zmax))
    top_clip = bounding_box_zmax
  
  # disturb values a tiny bit to avoid singularities
  eps1 = 1e-10  # 1.234e-2
  eps2 = 1e-10  # 5.432e-2
  bottom_clip += eps1
  top_clip -= eps2

  z_samples = np.linspace(bottom_clip, top_clip, n_loops)

  # initialize list of empty lists for rings
  loops = []
  for i in range(n_loops):
    loops.append([])

  debug = False

  print("Create {} loops for z in [{},{}] from the mesh {}.".format(n_loops, bottom_clip, top_clip, input_filename))
  # loop over z samples
  for loop_no,z_value in enumerate(z_samples):
      
    print("Loop no {}/{}".format(loop_no,len(z_samples)))

    # compute all intersecting line segments that lie in the surface and have the specified z_value
    create_loop(z_value, stl_mesh, loops[loop_no])
    
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
      
    new_loop = order_loop(loop, first_point)
            
    # store new loop 
    loops[loop_no] = list(new_loop)
    
    # adjust z value to undo artificial distortion
    for i in range(len(loops[loop_no])):
      loops[loop_no][i][2] = z_samples[loop_no]
        
  if debug:
    print("----------------")
    print("")
    print("loops: ",loops)
          
  print("The following loops have been extracted:")
  for (loop,z_value) in zip(loops,z_samples):
    print("at z = {:0.3f} comprising {} segments".format(z_value,len(loop)))
          
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

def create_boundary_points(input_filename, bottom_clip, top_clip, n_loops, n_points):
  """ 
  This is a top-level function that performs all the steps to create the initial boundary points of the whole mesh. 
  It is called by the C++ implementation.
  :param input_filename: filenamem of either an STL mesh or a geomdl B-spline Surface stored as pickle
  :param bottom_clip:
  :param top_clip:
  """
  
  # if the file contains a spline curve, try to load it and call the function from spline_surface
  if ".pickle" in input_filename:
    print("Interpreting \"{}\" as pickle file containing a NURBS surface.".format(input_filename))
    try:
      # check if file exists
      if not os.path.exists(input_filename):
        print("File \"{}\" does not exist".format(input_filename))
        
      f = open(input_filename,"rb")
      surface = pickle.load(f)
      return spline_surface.create_boundary_points(surface, bottom_clip, top_clip, n_loops, n_points)
    except:
      print("Error! Could not create boundary points from file \"{}\".".format(input_filename))
      quit()
  
  # algorithm for STL mesh
  print("Interpreting \"{}\" as STL file containing the surface.".format(input_filename))
  loops = create_rings(input_filename, bottom_clip, top_clip, n_loops, False)   # last argument is if debugging output should be written, set to False
  boundary_points, lengths = stl_create_mesh.rings_to_boundary_points(loops, n_points)
  boundary_points = stl_create_mesh.boundary_point_loops_to_list(boundary_points)
  
  return boundary_points

def create_point_marker(point, markers, size=None):
  if size is None:
    size = 0.02
  diag0 = np.array([-size,-size,-size])
  diag1 = np.array([size,-size,-size])
  diag2 = np.array([-size,size,-size])
  diag3 = np.array([size,size,-size])
  diag4 = np.array([-size,-size,size])
  diag5 = np.array([size,-size,size])
  diag6 = np.array([-size,size,size])
  diag7 = np.array([size,size,size])
  markers += [
    [point+diag0,point+diag3,point+diag1],[point+diag0,point+diag2,point+diag3],  # bottom
    [point+diag4,point+diag5,point+diag7],[point+diag4,point+diag7,point+diag6],  # top
    [point+diag0,point+diag1,point+diag5],[point+diag0,point+diag5,point+diag4],  # front
    [point+diag2,point+diag7,point+diag3],[point+diag2,point+diag6,point+diag7],  # back
    [point+diag2,point+diag0,point+diag4],[point+diag2,point+diag4,point+diag6],  # left
    [point+diag1,point+diag3,point+diag7],[point+diag1,point+diag7,point+diag5]  # right
  ]

def get_stl_mesh(input_filename):
  """
  Get an stl mesh object from the stl file given in input_filename, for use with function create_ring_section_mesh
  """
  return mesh.Mesh.from_file(input_filename)
  
def create_ring_section(input_filename, start_point, end_point, z_value, n_points):
  """
  Create a curve on the intersection of a horizontal plane given by z_value and the surface from the stl file.
  From nearest point to start_point to nearest point to end_point, the direction is such that the length of the curve is minimal (there are 2 possible orientations cw/ccw)
  :param input_filename: file name of an stl file that contains the closed surface mesh of the muscle, aligned with the z-axis, or alternatively a pickle file of the surface
  :param start_point: the line starts at the point on the surface with given z_value, that is the nearest to start_point
  :param end_point: the line ends at the point on the surface with given z_value, that is the nearest to end_point
  :param z_value: the z level of the line on the surface
  :param n_points: number of points on the boundary
  :return: list of points
  """
  
  # if the file contains a NURBS surface, try to load it and call the function from spline_surface
  if ".pickle" in input_filename:
    try:
      f = open(input_filename,"rb")
      surface = pickle.load(f)
    except:
      print("Error! Could not open file \"{}\" for reading spline surface.".format(input_filename))

    debugging_points = []
    result = spline_surface.create_ring_section(surface, start_point, end_point, z_value, n_points, debugging_points)
    
    # debugging output
    if False:
      level = 0
      rank_no = z_value
      filename = "00_{}_pass".format(start_point)
      if len(debugging_points) != 0:
        stl_debug_output.output_points(filename, rank_no, level, debugging_points, 0.02)

      filename = "00_{}_start_end".format(start_point)
      stl_debug_output.output_points(filename, rank_no, level, [start_point, end_point], 0.1)

      filename = "00_{}_points".format(start_point)
      stl_debug_output.output_points(filename, rank_no, level, result, 0.05)

    return result
    
  # else interpret the file as stl mesh and use the stl mesh algorithm
  print("Interpreting \"{}\" as STL file containing the surface.".format(input_filename))
  stl_mesh = get_stl_mesh(input_filename)
  return create_ring_section_mesh(stl_mesh, start_point, end_point, z_value, n_points)
      
def create_ring_section_mesh(stl_mesh, start_point, end_point, z_value, n_points):
  """
  Create a curve on the intersection of a horizontal plane given by z_value and the surface from the stl file.
  From nearest point to start_point to nearest point to end_point, the direction is such that the length of the curve is minimal (there are 2 possible orientations cw/ccw)
  :param stl_mesh: stl mesh that contains the closed surface mesh of the muscle, aligned with the z-axis, can be retrieved by get_stl_mesh
  :param start_point: the line starts at the point on the surface with given z_value, that is the nearest to start_point
  :param end_point: the line ends at the point on the surface with given z_value, that is the nearest to end_point
  :param z_value: the z level of the line on the surface
  :param n_points: number of points on the boundary
  :return: list of points
  """
  
  print("create_ring_section, start_point={}, end_point={}, z_value={}, n_points={}".format(start_point, end_point, z_value, n_points))
  
  debug = False                 # set this to true to enable debugging output
  write_output_mesh = False    # set this to true to output 4 stl meshes that explain the algorithm
  
  # create a full loop at the given z_value
  loop = []
  create_loop(z_value, stl_mesh, loop)
  
  new_loop = order_loop(loop, loop[0][0])
  
  if debug:
    print("new_loop: {}".format(new_loop))
  
  # determine edge which is closest to start_point
  min_distance_start_point = None
  min_distance_end_point = None
  
  near_start_point = np.array(start_point)
  near_end_point = np.array(end_point)
  loop_start_point = None
  loop_end_point = None
  edge_index_end = 0
  edge_index_start = 0
  
  # loop over edges
  for i,edge in enumerate([[new_loop[i], new_loop[(i+1)%len(new_loop)]] for i in range(len(new_loop))]):
    p0 = np.array(edge[0])
    p1 = np.array(edge[1])
    u = -p0 + p1   # edge
    
    if debug:
      print("edge {}, u: {}".format(edge, u))
    
    # determine loop_start_point as the point in the loop that is closest to the given near_start_point
    if u.dot(u) < 1e-10:
      t_start = -1
    else:
      t_start = (near_start_point - p0).dot(u) / u.dot(u)

    if debug:
      print("t_start: {}".format(t_start))
      
    if t_start >= 0 and t_start <= 1:
      plumb_foot_point_start = p0 + t_start * u
      
      if debug:
        print("")
        print("start point found, plumb_foot_point_start: ",plumb_foot_point_start)
        print("")
      
      distance_start_point = np.linalg.norm(near_start_point - plumb_foot_point_start)
      if min_distance_start_point is None or distance_start_point < min_distance_start_point:
        min_distance_start_point = distance_start_point
        loop_start_point = plumb_foot_point_start
        edge_index_start = (i+1)%len(new_loop)
        
    if min_distance_start_point is None or np.linalg.norm(near_start_point - p0) < min_distance_start_point:
      min_distance_start_point = np.linalg.norm(near_start_point - p0)
      loop_start_point = p0
      edge_index_start = i
      
      if debug:
        print("take p0 ({}), new min_distance_start_point: {}".format(p0,min_distance_start_point))
      
    if min_distance_start_point is None or np.linalg.norm(near_start_point - p1) < min_distance_start_point:
      min_distance_start_point = np.linalg.norm(near_start_point - p1)
      loop_start_point = p1
      edge_index_start = (i+1)%len(new_loop)

      if debug:
        print("take p1 ({}), new min_distance_start_point: {}".format(p1,min_distance_start_point))
        
    # determine loop_end_point as the point in the loop that is closest to the given near_end_point
    if u.dot(u) < 1e-10:
      t_end = -1
    else:
      t_end = (near_end_point - p0).dot(u) / u.dot(u)
    
    if debug:
      print("t_end: {}".format(t_end))
    
    if t_end >= 0 and t_end <= 1:
      plumb_foot_point_end = p0 + t_end * u
        
      if debug:
        print("")
        print("end point found, plumb_foot_point_end: ",plumb_foot_point_end)
        print("")
      
      distance_end_point = np.linalg.norm(near_end_point - plumb_foot_point_end)
      if min_distance_end_point is None or distance_end_point < min_distance_end_point:
        min_distance_end_point = distance_end_point
        loop_end_point = plumb_foot_point_end
        edge_index_end = (i+1)%len(new_loop)
      
    if min_distance_end_point is None or np.linalg.norm(near_end_point - p0) < min_distance_end_point:
      min_distance_end_point = np.linalg.norm(near_end_point - p0)
      loop_end_point = p0
      edge_index_end = i
      
      if debug:
        print("take p0 ({}), new min_distance_end_point: {}".format(p0,min_distance_end_point))
      
    if min_distance_end_point is None or np.linalg.norm(near_end_point - p1) < min_distance_end_point:
      min_distance_end_point = np.linalg.norm(near_end_point - p1)
      loop_end_point = p1
      edge_index_end = (i+1)%len(new_loop)
      
      if debug:
        print("take p0 ({}), new min_distance_end_point: {}".format(p0,min_distance_end_point))

  # here edge_index_start and edge_index_end are the indices of the first loop points after the start point (loop_start_point) and end point (loop_end_point), which usually lie on an edge and not at a node of the loop
  
  if debug:
    print("start_point: {}".format(start_point))
    print("end_point: {}".format(end_point))
    print("min_distance_start_point: {}".format(min_distance_start_point))
    print("min_distance_end_point: {}".format(min_distance_end_point))
    print("loop_start_point: {}".format(loop_start_point))
    print("loop_end_point: {}".format(loop_end_point))
  
  # determine length of loop between loop_start_point and loop_end_point, from start to end (distance0) and from end to start (distance1), always in forward direction 
  # iterate over points of loop
  previous_point = np.array(new_loop[0])
  point_index = 1
  first_section = 0   # 0 = not yet started, 1 = currently running, 2 = finished
  second_section = 0   # 0 = not yet started, 1 = currently running, 2 = finished
  distance0 = 0
  distance1 = 0
  start_point_index = edge_index_start
  end_point_index = edge_index_end
  
  while True:
    current_point = np.array(new_loop[point_index])
    # always consider the edge between previous_point and current_point
    
    if debug:
      print("previous_point: {}, current_point: {}, first: {}, second: {}, d: {},{}".format(previous_point, current_point, first_section, second_section, distance0, distance1))
    
    if second_section == 1:
      distance1 += np.linalg.norm(current_point - previous_point)
      
      if debug:
        print("  add to distance1 {}, new: {}".format(np.linalg.norm(current_point - previous_point), distance1))
      
    # if loop_end_point lies on the line between previous_point and current_point
    if point_index == edge_index_end:
      
      if debug:
        print("line with loop_end_point found at point_index {}".format(point_index))
      
      if first_section == 1:
        distance0 += np.linalg.norm(loop_end_point - previous_point)
        
        if debug:
          print("add to distance0 {}, new: {}".format(np.linalg.norm(loop_end_point - previous_point), distance0))
        first_section = 2
      
      if second_section == 0:
        second_section = 1
        distance1 = np.linalg.norm(current_point - loop_end_point)
              
        if debug:
          print("set distance1: {}".format(distance1))
    
    if first_section == 1:
      distance0 += np.linalg.norm(current_point - previous_point)
      
      if debug:
        print("  add to distance0 {}, new: {}".format(np.linalg.norm(current_point - previous_point), distance0))
      
    # if loop_start_point lies on the line between previous_point and current_point
    if point_index == edge_index_start:
      
      if debug:
        print("line with loop_start_point found at point_index {}".format(point_index))
      
      if first_section == 0:
        first_section = 1
        distance0 = np.linalg.norm(current_point - loop_start_point)
        
        if debug:
          print("set distance0: {}".format(distance0))
        
      if second_section == 1:
        distance1 -= np.linalg.norm(current_point - previous_point)
        distance1 += np.linalg.norm(loop_start_point - previous_point)
        
        if debug:
          print("add to distance1: {}, new: {}".format(np.linalg.norm(loop_start_point - previous_point)-np.linalg.norm(current_point - previous_point), distance1))
        
        second_section = 2
      
    if first_section == 2 and second_section == 2:
      break
      
    # go to next point
    point_index = (point_index+1) % len(new_loop)
    previous_point = current_point
  
    if debug:
      print("distance0: {}, distance1: {}".format(distance0, distance1))
            
  direction_reversed = False
  # here we have the distance between loop_start_point and loop_end_point following the ring forward (distance0) and backward (distance1)
  if distance0 > distance1:
    
    direction_reversed = True
    if debug:
      print("reverse direction")
    distance0,distance1 = distance1,distance0
    loop_start_point,loop_end_point = loop_end_point,loop_start_point
    start_point_index,end_point_index = end_point_index,start_point_index

  if debug:  
    print("loop_start_point: {}, loop_end_point: {}".format(loop_start_point, loop_end_point))
    print("start_point_index: {}, end_point_index: {}, length of new_loop: {}".format(start_point_index, end_point_index, len(new_loop)))
    print("distance0: {}".format(distance0))
  
  element_length = (float)(distance0) / (n_points-1)
  
  if debug:
    print("{} points, element_length: {}".format(n_points, element_length))
  
  markers_result = []
  # sample points in ring section
  point_index = start_point_index
  
  result_points = [loop_start_point]
  
  previous_point = loop_start_point
  current_point = np.array(new_loop[point_index])
  
  t = 0.0
  while True:
    
    distance_current_previous = np.linalg.norm(current_point - previous_point)
    while t + distance_current_previous > element_length:
      alpha = (element_length-t) / distance_current_previous
      new_point = (1.-alpha)*previous_point + alpha*current_point
      result_points.append(new_point)
      
      markers_result.append([previous_point, new_point, current_point])
      
      previous_point = new_point
      t = 0.0
      distance_current_previous = np.linalg.norm(current_point - previous_point)
      
    t += distance_current_previous  # < element_length
    
    # if end is reached
    if point_index == end_point_index:
      break
      
    # go to next point
    previous_point = current_point
    point_index = (point_index+1) % len(new_loop)
    current_point = np.array(new_loop[point_index])
  
    # if end is reached
    if point_index == end_point_index:
      
      if debug:
        print("reached end point")
      current_point = loop_end_point
  
  if not np.allclose(result_points[-1], loop_end_point):    
    if debug:
      print("add end point")
    result_points.append(loop_end_point)
  
    
  # result_points = [loop_start_point]
  # element_distance = 0
  # previous_element_distance = 0
    
  # previous_point = loop_start_point
  # current_point = np.array(new_loop[point_index])
    
  # element_distance += np.linalg.norm(current_point - previous_point)
  # while True:
    # print("element_distance: {}".format(element_distance))
    # #current_point = np.array(new_loop[point_index])
      
    # #element_distance += np.linalg.norm(current_point - previous_point)
      
    # if element_distance > element_length:
      # t = (element_length-previous_element_distance) / (element_distance-previous_element_distance)
      # new_point = previous_point + t*(current_point - previous_point)
      # result_points.append(new_point)
      # print("    point {} between {} and {}, t={}".format(new_point, previous_point, current_point, t))
      # previous_point = new_point
      # element_distance = 0
      # distance_current_previous = np.linalg.norm(current_point - previous_point)
      # previous_element_distance = element_distance
      # element_distance += distance_current_previous
    # else:
        
      # if point_index == end_point_index:
        # break
      # # go to next point
      # point_index = (point_index+1) % len(new_loop)
      # print("    go to next point index {}, end_point_index: {}".format(point_index, end_point_index))
      # previous_point = current_point
      # current_point = np.array(new_loop[point_index])
      # distance_current_previous = np.linalg.norm(current_point - previous_point)
      # previous_element_distance = element_distance
      # element_distance += distance_current_previous
    
  if debug:
    print ("result contains {} points".format(len(result_points)))
      
  # transform points from numpy array to simple list
  result = []
  if direction_reversed:
    for point in reversed(result_points):
      result.append(point.tolist())
  else:
    for point in result_points:
      result.append(point.tolist())
        
  if write_output_mesh:

    # create markers
    markers_loop = []
    for point in new_loop:
      create_point_marker(point, markers_loop)
      
    markers_start_end = []
    markers_loop_start_end = []
    create_point_marker(start_point, markers_start_end, 0.1)
    create_point_marker(end_point, markers_start_end, 0.1)
    create_point_marker(loop_start_point, markers_loop_start_end, 0.08)
    create_point_marker(loop_end_point, markers_loop_start_end, 0.16)
        
    
    for point in result_points:
      create_point_marker(point, markers_result, 0.01)
      
    #---------------------------------------
    # Create the mesh
    out_mesh = mesh.Mesh(np.zeros(len(markers_start_end), dtype=mesh.Mesh.dtype))
    for i, f in enumerate(markers_start_end):
      out_mesh.vectors[i] = f
      #for j in range(3):
        #print("set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
        
    #out_mesh.update_normals()

    outfile = "mesh_start_end.stl"
    out_mesh.save(outfile)
    print("saved {} triangles to \"{}\" (start,end)".format(len(markers_start_end),outfile))

    #---------------------------------------
    # Create the mesh
    out_mesh = mesh.Mesh(np.zeros(len(markers_loop_start_end), dtype=mesh.Mesh.dtype))
    for i, f in enumerate(markers_loop_start_end):
      out_mesh.vectors[i] = f
      #for j in range(3):
        #print("set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
        
    #out_mesh.update_normals()

    outfile = "mesh_loop_start_end.stl"
    out_mesh.save(outfile)
    print("saved {} triangles to \"{}\" (loop start,end)".format(len(markers_loop_start_end),outfile))

    #---------------------------------------
    out_mesh = mesh.Mesh(np.zeros(len(markers_loop), dtype=mesh.Mesh.dtype))
    for i, f in enumerate(markers_loop):
      out_mesh.vectors[i] = f
      #for j in range(3):
        #print("set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
        
    #out_mesh.update_normals()

    outfile = "mesh_loop.stl"
    out_mesh.save(outfile)
    print("saved {} triangles to \"{}\" (markers_loop)".format(len(markers_loop),outfile))

    #---------------------------------------
    out_mesh = mesh.Mesh(np.zeros(len(markers_result), dtype=mesh.Mesh.dtype))
    for i, f in enumerate(markers_result):
      out_mesh.vectors[i] = f
      #for j in range(3):
        #print("set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
        
    #out_mesh.update_normals()

    outfile = "mesh_result.stl"
    out_mesh.save(outfile)
    print("saved {} triangles to \"{}\" (start,end)".format(len(markers_result),outfile))

  if debug:
    print("result:{}".format(result))
  return result
    
