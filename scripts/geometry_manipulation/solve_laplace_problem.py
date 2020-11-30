#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file gets imported by stl_create_mesh.py

import sys, os
import numpy as np

import scipy
import scipy.spatial
import scipy.integrate
import scipy.optimize

def solve_laplace_problem(n_points, points, n_original_points, original_points, original_point_phi_value, triangulation_type, parametric_space_shape, point_indices_list, triangle_list, debugging_stl_output, stl_triangle_lists):
  """
  This is a helper function for stl_create_mesh.create_planar_mesh, defined in create_planar_mesh.py
  Solve the Laplace problem / potential flow to construct the harmonic map from world domain to parameter domain
  
  :param n_points: number of points on the boundary of the resulting mesh
  :param points: a vector of 3D points with shape (n_points x 3), each row is one point
  :param n_original_points: size of original_points
  :param original_points: a copy of points made before additional boundary points were added
  :param parametric_space_shape: 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
  :param point_indices_list: a list of the indices into the points array for each triangle of the triangulation
  :param triangle_list: the resulting triangles with their points
  :param debugging_stl_output: if list should be filled with STL triangles that can be output to a STL mesh for debugging
  :param stl_triangle_lists: the debugging lists: [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,markers_grid_points_parametric_space, markers_grid_points_world_space]
  :return: u,v the solution vectors for the two harmonic maps
  """
  
  debug = False   # enable debugging output
  if debugging_stl_output:
    [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,\
      markers_grid_points_parametric_space, markers_grid_points_world_space] = stl_triangle_lists

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
    if parametric_space_shape == 0 or parametric_space_shape == 3:  # unit circle
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
      
  return u,v
  
