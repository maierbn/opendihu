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

def create_slice_tringulation(triangulation_type, modify_phi, n_points, points, max_area_factor):
  """
  This is a helper function for stl_create_mesh.create_planar_mesh, defined in create_planar_mesh.py
  Create a 2d triangulation on a slice using different methods.
  :param triangulation_type: 0 = scipy, 1 = triangle, 2 = center pie (2 is best), 3 = minimized distance
  :param modify_phi: if the grid should be rotated
  :param n_points: number of points on the boundary of the resulting mesh
  :param points: a vector of 3D points with shape (n_points x 3), each row is one point
  :param max_area_factor: only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
  :return: point_indices_list,triangle_list,n_points
    point_indices_list: (output) a list of the indices into the points array for each triangle of the triangulation
    triangle_list: (output) the resulting triangles with their points
  """
  
  # extract z value 
  center_point = np.sum(points,axis=0)/len(points)
  z_value = center_point[2]
  
  # get information about extent of points
  (max_x,max_y,max_z) = np.max(points,axis=0)
  (min_x,min_y,min_z) = np.min(points,axis=0)
  extent_x = max_x - min_x
  extent_y = max_y - min_y
  
  # store points, because later they will be overwritten by adding new points from the triangulation
  original_points = np.array(points)
  n_original_points = n_points
  n_points_per_face = (int)(n_original_points/4)
  n_regular_grid_boundary_points = n_original_points
  
  # project points on xy=z_value plane
  projected_points = []
  for point in points:
    projected_points.append(np.array([point[0], point[1]]))
  
  projected_points = np.reshape(projected_points, (-1,2))

  # define helper variables with default values that are only later used when triangulation_type == 2
  def get_modified_phi(phi_in): 
    phi_out = phi_in
    return phi_out
  determine_additional_points_on_ring = False
  n_additional_points_on_ring = 0
  original_point_phi_value = []
  
  debug = False   # enable debugging output
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
      
    triangulated_projected_points = np.array(triangulation['vertices'])
    
    # transform projected points back to 3D points
    points = []
    for projected_point in triangulated_projected_points:
      points.append(np.array([projected_point[0], projected_point[1], z_value]))
    
    # update n_points
    n_points = len(points)
    points = np.reshape(points, (-1,3))
    
    # determine the phi angle in the circle of the current point
    
    for original_point_no in range(n_original_points):
      boundary_point = points[original_point_no]
    
      phi = float(original_point_no) / n_original_points * 2 * np.pi
      original_point_phi_value.append(phi)
    
    # add additional points on ring
    # settings
    determine_additional_points_on_ring = True
    rescale_phi = True
    
    # normal implementation without searching for additional boundary points on ring that the triangulation created
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
        # check if this point lies between two boundary points
        # loop over initial boundary points
        for boundary_point_index in range(n_original_points):
          boundary_point1 = points[boundary_point_index]
          boundary_point2 = points[(boundary_point_index+1)%n_original_points]
          
          v1 = -point + boundary_point1
          v2 = -point + boundary_point2
          v3 = -boundary_point1 + boundary_point2
          
          v1 = v1 / np.linalg.norm(v1)
          v2 = v2 / np.linalg.norm(v2)
          
          # if the point lies on the line between the two boundary points
          if abs(np.linalg.norm(np.cross(v1,v2))) < 1e-3:
            if abs(v3[0]) < abs(v3[1]):            
              alpha = (point[1] - boundary_point1[1]) / v3[1]
            else:
              alpha = (point[0] - boundary_point1[0]) / v3[0]
            
            if alpha > 1.0 or alpha < 0.0:
              #print("alpha: {} continue".format(alpha))
              continue
            
            phi = original_point_phi_value[boundary_point_index] + alpha * \
               (original_point_phi_value[(boundary_point_index+1)%n_original_points] - original_point_phi_value[boundary_point_index])
            original_point_phi_value.append(phi)
        
            #print("{} alpha: {}, phi: {} {} -> {}".format(point_index, alpha, original_point_phi_value[boundary_point_index], original_point_phi_value[(boundary_point_index+1)%n_original_points], phi))
        
            new_points.append(point)
            additional_points_on_ring.append(point_index)
            current_point_is_on_boundary = True
            break
        
        if not current_point_is_on_boundary:
          interior_points.append(point)
      
      # store points such that points = [<previous original points>, <newly determined points on the boundary>, <rest of points>]
      
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
            
            # count number of boundary points before old point_index
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
      
      # points has the following structure: [<list of original boundary points>, <list of new boundary points>, <list of interior points>]
      # original_points has the following structure: [<list of original boundary points>, <list of new boundary points>]
      points = np.array(new_points)
      triangle_list = points[point_indices_list]
      n_regular_grid_boundary_points = n_original_points
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
    
      #print("previous boundary points: {}, n_original_points: {}, n_additional_points_on_ring: {}, interior: {}, n_points: {}".\
      #  format(n_original_points-n_additional_points_on_ring, n_original_points, n_additional_points_on_ring, n_points-n_original_points, n_points))
      #print("additional_points_on_ring: {}".format(additional_points_on_ring))
      
      # setup map between parameter space regular grid in the circle and the transformed parameter space grid with the additional boundary points
      # this is done by defining a map for phi
      # map from phi to phi
      
      def get_modified_phi(phi_in):
        
        # normalize phi to [0,2*np.pi)
        if phi_in < 0:
          phi_in += 2*np.pi
        
        # determine position of phi between regular grid boundary points
        phi_increment = (2*np.pi) / n_regular_grid_boundary_points
        previous_boundary_point_index = (int)(phi_in / phi_increment)
        
        # determine factor between previous and next boundary point
        alpha = (phi_in - previous_boundary_point_index*phi_increment) / phi_increment
        
        # determine positions of phi in the new boundary points
        next_phi_value = 2*np.pi
        if previous_boundary_point_index+1 < len(original_point_phi_value):
          next_phi_value = original_point_phi_value[previous_boundary_point_index+1]
          
        previous_phi_value = original_point_phi_value[previous_boundary_point_index]
        
        # compute phi value with alpha between new boundary points
        phi_out = previous_phi_value + alpha * (next_phi_value - previous_phi_value)
        
        #print("phi_in: {}, phi_increment: {}, previous_boundary_point_index:{} [{},{}], alpha:{} new:[{},{}], phi_out: {}".format(phi_in, phi_increment, previous_boundary_point_index, previous_boundary_point_index*phi_increment, (previous_boundary_point_index+1)*phi_increment, alpha,\
        #  previous_phi_value, next_phi_value, phi_out))
        
        return phi_out
    
  elif triangulation_type == 2 or triangulation_type == 3:
    # 2: simple custom triangulation with triangles around one center point in CoG
    # 3: custom triangulation with triangles around point for which distance is minimized

    # compute the center point by minimizing the distances to the boundary points
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
      
      # compute the rotation angle when iterating over all connection vectors between center and boundary point
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
    
  #print("  number of projected points: ",len(projected_points),", number of initial triangles: ", len(point_indices_list))
  return point_indices_list, triangle_list, n_points, points, n_original_points, original_points, original_point_phi_value, get_modified_phi, n_regular_grid_boundary_points, extent_x, extent_y, n_additional_points_on_ring, determine_additional_points_on_ring
  
