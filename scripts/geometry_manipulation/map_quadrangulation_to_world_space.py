#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file gets imported by stl_create_mesh.py

import sys, os
import numpy as np

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
  
def transform_to_world_space(x, y, triangles_parametric_space, triangle_list, parametric_space_shape):
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
      if parametric_space_shape == 0 or parametric_space_shape == 3:  # unit circle
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
  
def map_quadrangulation_to_world_space(u, v, n_grid_points_x, n_grid_points_y, parametric_space_shape, modify_phi, original_point_phi_value, get_modified_phi, n_regular_grid_boundary_points, point_indices_list, triangle_list, extent_x, points, debugging_stl_output, stl_triangle_lists):
  """
  This is a helper function for stl_create_mesh.create_planar_mesh, defined in create_planar_mesh.py
  Map a quadrangulation on the parameter domain to the world domain (2D slice)
  :param u: harmonic map for x direction (solution of the Laplace problem)
  :param v: harmonic map for x direction (solution of the Laplace problem)
  :param n_grid_points_x: number of grid points in x direction of the final quadrangulation
  :param n_grid_points_y: number of grid points in y direction of the final quadrangulation
  :param parametric_space_shape: 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
  
  :param point_indices_list: a list of the indices into the points array for each triangle of the triangulation
  :param triangle_list: the resulting triangles with their points
  :param debugging_stl_output: if list should be filled with STL triangles that can be output to a STL mesh for debugging
  :param stl_triangle_lists: the debugging lists: [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,markers_grid_points_parametric_space, markers_grid_points_world_space]
  :return: grid_points_world_space, grid_points_parametric_space, the mapped grid points of the quadrangulation
  """
  
  debug = False   # enable debugging output
  if debugging_stl_output:
    [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,\
      markers_grid_points_parametric_space, markers_grid_points_world_space] = stl_triangle_lists

  
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
    center_point = np.sum(points,axis=0)/len(points)
    z_value = center_point[2]
  
    x_offset = center_point[0] + extent_x*1.5
    y_offset = center_point[1]
    scale = 10.0
    for triangle_parametric_space in triangles_parametric_space:
      out_triangle_parametric_space = []
      for p in triangle_parametric_space:
        out_triangle_parametric_space.append(np.array([p[0,0]*scale+x_offset, p[1,0]*scale+y_offset, z_value]))
        
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
      # if there were additional boundary points from the triangulation that also got Dirichlet BC, and modifiy_phi is set to true, compute world points differently
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
    
  #distances_current_loop,relative_distances_current_loop = compute_mean_distances(grid_points_world_space)
  #print("transformed, std: {}".format(np.std(relative_distances_current_loop))
  
  return grid_points_world_space, grid_points_parametric_space, grid_points_parametric_space_modified
