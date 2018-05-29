#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script creates a hexaeder mesh out of the rings given by create_rings.py. 
# It reads in a pickle file `rings` and creates a pickle file `mesh` containing 
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
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import struct
import stl
from stl import mesh
import pickle
import scipy
import scipy.spatial
import scipy.linalg
import scipy.integrate

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
    print "xi: {}, {}".format(xi1,xi2)
  
  # point is inside triangle iff conditions for xi1,xi2 are fulfilled
  tol = 1e-14
  condition = (-tol <= xi1 <= 1+tol and -tol <= xi2 <= 1+tol and xi1 + xi2 <= 1+tol)
  
  if debug:
    if condition:
      print "triangle: ({},{}) ({},{}) ({},{}), point: ({},{})".format(point0[0], point0[1], point1[0], point1[1], point2[0], point2[1], point[0], point[1])
    
      print "dudxi: \n{}".format(dudxi)
      print "dxidu: \n{}".format(dxidu)
      print "p-p0: \n{}".format(point-point0)
      
      print "point found"
      print ""
  
  return (condition, xi)
  
# constant parameters
triangulation_type = 1  # 0 = scipy, 1 = triangle, 2 = custom (1 is best)
parametric_space_shape = 3   # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
max_area_factor = 2.    # only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
show_plot = False

n_points_x = 4    # number of points per dimension for the initial triangulation, assuming a cube reference area (4*n_points_x on circumference for unit circle)

n_grid_points_x = n_points_x+1   # grid width of generated 2d mesh
n_grid_points_y = n_points_x+1

if parametric_space_shape == 0:  # for unit circle 
  n_grid_points_y = 20
  
  
if len(sys.argv) < 2:
  print "usage: ./create_mesh.py [<triangulation_type> [<parametric_space_shape> [<n_points_x> [<n_grid_points_x>]]]]"
  sys.exit(0)

if len(sys.argv) >= 2:
  triangulation_type = int(sys.argv[1])
  
if len(sys.argv) >= 3:
  parametric_space_shape = int(sys.argv[2])
  
if len(sys.argv) >= 4:
  n_points_x = int(sys.argv[3])    
  n_grid_points_x = n_points_x+1   # grid width of generated 2d mesh
  n_grid_points_y = n_points_x+1
  
if len(sys.argv) >= 5:
  n_grid_points_x = int(sys.argv[4])
  n_grid_points_y = n_grid_points_x
  
print "triangulation_type: {}".format(triangulation_type)
print "parametric_space_shape: {}".format(parametric_space_shape)
print "n_points_x: {}".format(n_points_x)
print "n_grid_points_x: {}".format(n_grid_points_x)
print "n_grid_points_y: {}".format(n_grid_points_y)
  
# read in loops
with open('rings', 'rb') as f:
  loops = pickle.load(f)
  
debug = False  

print "{} loops".format(len(loops))

# iterate over loops of points
sorted_loops = []
lengths = []
for loop in loops:

  if not loop:
    continue

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
    
  if debug:
    print "loop"
    print loop
    print ""
    
    print "sorted_loop"
    print sorted_loop
    print ""
    print ""
    
  sorted_loops.append(sorted_loop)

# note: the first and last point in the loop is the same, because it represents the edges between the points

n_points = 4*n_points_x
  
# triangle lists for debugging output to stl files
out_triangulation_world_space = []
markers_border_points_world_space = []
out_triangulation_parametric_space = []
grid_triangles_world_space = []
grid_triangles_parametric_space = []
markers_grid_points_parametric_space = []
markers_grid_points_world_space = []
 
loop_grid_points = []  # list of grid point, for every slice, only contains loops that are not empty
  
debug = False
  
# sample loop with 4*n_points_x equidistant points
for loop_no,(loop,length) in enumerate(zip(sorted_loops,lengths)):
  
  if not loop:
    continue
    
  # for debugging only consider loop_no
  #if loop_no != 2:
  #  continue
    
  print ""
  print "Loop {}/{} with {} points, length: {}".format(loop_no, len(sorted_loops), len(loop), length)
  
  border_points = []
  n_points = 4*n_points_x
  h = float(length) / n_points
  
  # loop over points of loop
  previous_loop_point = None
  for loop_point in loop:
    
    # for the start of the loop
    if previous_loop_point is None:
      previous_loop_point = loop_point
      
      # add starting point of loop as first border points
      border_points.append(loop_point)
      if debug:
        print "add start point: ",border_points
      
      t_previous_loop_point = 0
      t_next_border_point = h
      continue
    
    if debug:
      print "current previous_loop_point: {}, t_previous_loop_point: {}".format(previous_loop_point,t_previous_loop_point) 
    
    # compute current edge
    edge = -previous_loop_point + loop_point 
    edge_length = np.linalg.norm(edge) 
    
    n_on_edge = 0
    
    # collect all border points on current edge
    while t_next_border_point <= t_previous_loop_point + edge_length:
      
      border_point = previous_loop_point + edge * (t_next_border_point - t_previous_loop_point) / edge_length
      border_points.append(border_point)
      
      if debug:
        print "border point {}, t previous: {}, edge length: {}, next: {}".format(border_point, t_previous_loop_point, edge_length, t_next_border_point)
      
      t_next_border_point += h
      if debug:
        n_on_edge += 1
    
    if debug:
      print "n_on_edge: ",n_on_edge
      
      
    # move on to next edge
    t_previous_loop_point += edge_length
    previous_loop_point = loop_point
  
  # handle edge back to first point
  loop_point = loop[0]

  # compute current edge
  edge = -previous_loop_point + loop_point 
  edge_length = np.linalg.norm(edge) 
  
  # collect all border points on current edge
  while t_previous_loop_point + edge_length >= t_next_border_point:
    
    border_point = previous_loop_point + edge * (t_next_border_point - t_previous_loop_point)
    border_points.append(border_point)
    
    if debug:
      print "border point {}, t previous: {}, edge length: {}, next: {}".format(border_point, t_previous_loop_point, edge_length, t_next_border_point)
    
    t_next_border_point += h
    
  
  # if there were too many points collected, due to rounding errors
  if len(border_points) > n_points:
    if debug:
      print "too many points: {}, n_points: {}".format(len(border_points), n_points)
    border_points = border_points[:n_points]

  
  # triangulate surface in world space
  if debug:
    print "border points: ",len(border_points)
    print border_points
  
  print "n border points: ",len(border_points)
  
  points = np.reshape(border_points,(n_points,3))
  
  if debug:  
    print ""
    print "points:"
    print points

  # for debugging create markers at border points
  for point in points:
    size = 0.2
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
  
  # project points on xy=z_value plane
  projected_points = []
  for point in points:
    projected_points.append(np.array([point[0], point[1]]))
  
  projected_points = np.reshape(projected_points, (-1,2))
  
  if debug:
    print ""
    print "projected_points:"
    print projected_points
  
  if triangulation_type == 0:
    # delaunay triangulation of scipy, does not add new points but is not constrained (works for convex domains)
  
    triangulation = scipy.spatial.Delaunay(projected_points)
    point_indices_list = triangulation.simplices
    triangle_list = points[point_indices_list]
    
  elif triangulation_type == 1:
    # delaunay triangulation of triangle package, adds new points, is constrained (works for concave domains)
    
    import triangle   # sudo easy_install triangle

    # create delaunay triangulation of points
    segments = np.reshape([[i,i+1] for i in range(n_points)], (n_points,2))
    segments[n_points-1] = np.array([n_points-1,0])
    
    data = {"vertices": projected_points, "segments": segments}
  
    max_area = extent_x * extent_y / max_area_factor
    if debug:
      print "maximum area: ",max_area
  
    #triangulation = triangle.triangulate(data, 'pq')
    triangulation = triangle.triangulate(data, 'pqa'+str(max_area))
    triangulated_projected_points = np.array(triangulation['vertices'])
    
    # transform projected points back to 3D points
    points = []
    for projected_point in triangulated_projected_points:
      points.append(np.array([projected_point[0], projected_point[1], z_value]))
    
    # update n_points
    n_points = len(points)
    points = np.reshape(points, (-1,3))
    
    point_indices_list = triangulation["triangles"]
    triangle_list = points[point_indices_list]
    
  elif triangulation_type == 2:
    # simple custom triangulation
    
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
    
    
    #print "point_indices_list:",point_indices_list
    #print "points:",points,points.shape
    triangle_list = []
    for point_index_triple in point_indices_list:
      triangle_list.append(points[point_index_triple])
    
    #triangle_list = points[point_indices_list]  ## doesn't work sometimes
    
  print "number of projected points: ",len(projected_points),", number of initial triangles: ", len(point_indices_list)
    
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
    print "n_dofs=",n_dofs,",n_original_points=",n_original_points
  
  # loop over boundary points, `original_points` are the points of the ring surface, `points` is a superset containing additional points created by the triangulation
  for original_point_no,original_point in enumerate(original_points):
    
    # get the position in reference space
    if parametric_space_shape == 0 or parametric_space_shape == 3:  # unit circle
      phi = float(original_point_no) / n_original_points * 2 * np.pi
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
      print "global stiffness matrix:"
      print global_stiffness

      print ""
      print "rhs vector (u):"
      print rhs_u

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
      print "rhs_non_dirichlet_u:",rhs_non_dirichlet_u
      print "global_stiffness_non_dirichlet:",global_stiffness_non_dirichlet
      print "non_dirichlet_u:",non_dirichlet_u
      print "non_dirichlet_v:",non_dirichlet_v
    
    u = np.concatenate([dirichlet_u[0:n_original_points], non_dirichlet_u])
    v = np.concatenate([dirichlet_v[0:n_original_points], non_dirichlet_v])
    
    if debug:
      print "u:",u
    
  # output solution
  for dof_no in range(n_dofs):
    
    # find the corresponding point in original_points
    point = points[dof_no]
    original_point_no = None
    for no,original_point in enumerate(original_points):
      if abs(point[0]-original_point[0]) < 1e-12 and abs(point[1]-original_point[1]) < 1e-12:
        original_point_no = no
        break
        
    if debug:
      print "dof {}, original point no: {}, dirichlet: ({},{}), solution: ({},{}), rhs: ({},{})".\
        format(dof_no, original_point_no, dirichlet_u[dof_no], dirichlet_v[dof_no], u[dof_no], v[dof_no], rhs_u[dof_no], rhs_v[dof_no])
      
  # store the triangles in parametric space
  triangles_parametric_space = []
  
  # loop over triangles indices
  for point_indices in point_indices_list:
    
    # loop over the dofs of the curren triangle and store (u,v) values of current triangle
    triangle_parametric_space = []
    for dof_no in point_indices:
      point_parametric_space = np.array([u[dof_no], v[dof_no]])
      
      triangle_parametric_space.append(point_parametric_space)
    
    # store triangle to list of triangles
    triangles_parametric_space.append(triangle_parametric_space)
    
    
  # create parametric space triangles for debugging output, move near ring and scale 10x
  x_offset = center_point[0] + extent_x*1.5
  y_offset = center_point[1]
  scale = 10.0
  for triangle_parametric_space in triangles_parametric_space:
    
    out_triangle_parametric_space = []
    for p in triangle_parametric_space:
      out_triangle_parametric_space.append(np.array([p[0]*scale+x_offset, p[1]*scale+y_offset, z_value]))
      
    out_triangulation_parametric_space.append(out_triangle_parametric_space)
    
  # now the mapping x -> u,v is computed
  # create new grid points on the ring that form a uniform mesh in parametric space

  n_grid_points = n_grid_points_x*n_grid_points_y
  grid_points_world_space = np.empty((n_grid_points,3))
  grid_points_parametric_space = np.empty((n_grid_points,2))
  
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
            if i >= j and i <= n_grid_points_x-1-j:   # bottom
          
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
        
      elif parametric_space_shape == 3:    # unit circle with adjusted grid points
        
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
              
        if n_grid_points_x%2 == 1 and i == int(n_grid_points_x/2) and j == int(n_grid_points_y/2):   # center point
          x = 0.
          y = 0.
          
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
          
      if xi_point is None:
        for n_tries in range(5):
          if parametric_space_shape == 0 or parametric_space_shape == 3:  # unit circle
            phi = np.arctan2(y,x)
            r = x / np.cos(phi)
            print " (x,y) = ({},{}), (phi,r)=({}deg,{})  (check: (x,y)=({},{}))".format(x,y,phi*180./np.pi,r,r*np.cos(phi),r*np.sin(phi))
            if abs(r) <= 1e-10:
              r_old = r
              r = 1e-4*np.sign(r)
              phi_old = phi
              phi += 1e-4
              x_old = x
              y_old = y
              x = r*np.cos(phi)
              y = r*np.sin(phi)
              print " [1] adjust grid point (x,y) = ({},{})->({},{}), phi={}->{}, r={}->{}".format(x_old,y_old,x,y,phi_old,phi,r_old,r)
            
            elif abs(r) >= 1.0-1e-10:
              r_old = r
              r = 0.99*np.sign(r)
              phi_old = phi
              phi += 1e-4
              x_old = x
              y_old = y
              x = r*np.cos(phi)
              y = r*np.sin(phi)
              print " [2] adjust grid point (x,y) = ({},{})->({},{}), phi={}->{}, r={}->{}".format(x_old,y_old,x,y,phi_old,phi,r_old,r)
            
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
              print " [3] adjust grid point (x,y) = ({},{})->({},{}), phi={}->{}, r={}->{}".format(x_old,y_old,x,y,phi_old,phi,r_old,r)
              
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
                print " xi = [{},{}] found after {} tries".format(xi[0],xi[1],n_tries+1)
                break
        
      if xi_point is None:
        
        print "Error: could not find triangle in parameter space for grid point (x,y) = ({},{}), r={}".format(x,y,r)
        print ""
        grid_points_world_space[j*n_grid_points_x+i] = np.array([0.0,0.0,0.0])
        
      if xi_point is not None:
          
        triangle_world_space = triangle_list[triangle_parameteric_space_no]
        p1 = triangle_world_space[0]
        p2 = triangle_world_space[1]
        p3 = triangle_world_space[2]
        
        point_world_space = (1 - xi[0] - xi[1])*p1 + xi[0]*p2 + xi[1]*p3
        grid_points_world_space[j*n_grid_points_x+i] = point_world_space
        grid_points_parametric_space[j*n_grid_points_x+i] = np.array([x,y])
      
  # store grid points in world space of current loop
  loop_grid_points.append(grid_points_world_space)
      
  # create triangles of new grid points mesh
  grid_point_indices_world_space = []
  
  # loop over grid points in parametric space
  for (j,y) in enumerate(np.linspace(0.0,1.0,n_grid_points_y)):
    phi = float(y) * (n_grid_points_y-1.0) / n_grid_points_y  * 2.*np.pi
    for (i,x) in enumerate(np.linspace(0.0,1.0,n_grid_points_x)):
  
      [x,y] = grid_points_parametric_space[j*n_grid_points_x+i]
  
      # for debugging create markers at grid points in parametric space
      point = np.array(np.array([x*scale+x_offset, y*scale+y_offset, z_value]))
      size = 0.2
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
      size = 0.2
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

      if parametric_space_shape == 1 or parametric_space_shape == 3:  # unit square 
        if i == n_grid_points_x-1 or j == n_grid_points_x-1:
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
      
      offset = np.array([x_offset, y_offset])
      p0 = np.concatenate([grid_points_parametric_space[j*n_grid_points_x+i]*scale+offset,np.array([z_value])])
      p1 = np.concatenate([grid_points_parametric_space[j*n_grid_points_x+(i+1)%n_grid_points_x]*scale+offset,np.array([z_value])])
      p2 = np.concatenate([grid_points_parametric_space[(j+1)%n_grid_points_y*n_grid_points_x+i]*scale+offset,np.array([z_value])])
      p3 = np.concatenate([grid_points_parametric_space[(j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x]*scale+offset,np.array([z_value])])
      
      grid_triangles_parametric_space.append([p0,p1,p3])
      grid_triangles_parametric_space.append([p0,p3,p2])

  # plot laplace solutions
  x = np.reshape(points[:,0], (-1))
  y = np.reshape(points[:,1], (-1))
  u_list = np.reshape(u, (-1))
  v_list = np.reshape(v, (-1))
  
  xw = np.reshape(grid_points_world_space[:,0], (-1))
  yw = np.reshape(grid_points_world_space[:,1], (-1))

  f, ax = plt.subplots(2,2)
  
  # u
  ax[0,0].tricontourf(x,y,u_list, 20) # 20 contour levels
  ax[0,0].triplot(x,y,point_indices_list,color='k')
  ax[0,0].plot(x,y, 'ko')
  ax[0,0].set_title('u, triangulation in world space')
  ax[0,0].set_aspect('equal')
  
  # v
  ax[0,1].tricontourf(x,y,v_list, 20) # 20 contour levels
  ax[0,1].triplot(x,y,point_indices_list,color='k')
  ax[0,1].plot(x,y, 'ko')
  ax[0,1].set_title('v, triangulation in world space')
  ax[0,1].set_aspect('equal')
  
  
  # parametric space
  ax[1,0].triplot(u_list,v_list,point_indices_list,color='k')
  ax[1,0].plot(u_list,v_list, 'ko')
  ax[1,0].set_title('triangulation in parametric space')
  ax[1,0].set_aspect('equal')
  
  # world space grid
  ax[1,1].triplot(xw,yw,grid_point_indices_world_space,color='k')
  ax[1,1].plot(xw,yw, 'ko')
  ax[1,1].set_title('new grid in world space')
  ax[1,1].set_aspect('equal')
  
  plt.savefig("out/harmonic_map_{}.png".format(loop_no))
  if show_plot:
    plt.show()
  plt.close()
  
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
  print "Grid per ring is {} x {}, {} rings: with this it is not possible to create quadratic elements (numbers must be odd)".format(n_grid_points_x, n_grid_points_y, n_loops)

# fill list of nodes
for grid_point_list in loop_grid_points:
  for node_position in list(grid_point_list):
    node_positions.append(list(node_position))
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
  
# create stl mesh of 3D mesh
out_3d_mesh_triangles = []
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
    seed_points.append(list(center_point))
  
# compute indices of bottom and top nodes. These are needed for setting boundary conditions
bottom_node_indices = range(n_grid_points_per_loop)
top_node_indices = range(n_grid_points_per_loop*(n_loops-1),n_grid_points_per_loop*n_loops)
  
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
with open('mesh', 'wb') as f:
  pickle.dump(data, f)
  
  
# write debugging output stl meshes
def write_stl(triangles, outfile, description):
  # create output mesh
  n_triangles = len(triangles)

  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
    #for j in range(3):
      #print "set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
      
  out_mesh.update_normals()

  out_mesh.save(outfile) #, mode=stl.Mode.ASCI
  print "saved {} triangles to \"{}\" ({})".format(n_triangles,outfile,description)

write_stl(markers_border_points_world_space,   "out/mesh_02_border_points_w.stl", "border points")
write_stl(out_triangulation_world_space,       "out/mesh_03_triangulation_w.stl", "triangulation world space")
write_stl(out_triangulation_parametric_space,  "out/mesh_04_triangulation_p.stl", "triangulation parametric space")
write_stl(grid_triangles_parametric_space,     "out/mesh_05_grid_triangles_p.stl","grid parametric space")
write_stl(markers_grid_points_parametric_space,"out/mesh_06_grid_points_p.stl",   "grid points parametric space")
write_stl(grid_triangles_world_space,          "out/mesh_07_grid_triangles_w.stl","grid world space")
write_stl(markers_grid_points_world_space,     "out/mesh_08_grid_points_w.stl",   "grid points world space")
write_stl(out_3d_mesh_triangles,               "out/mesh_09_3d_mesh_w.stl",       "3d mesh world space")



