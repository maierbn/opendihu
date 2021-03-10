#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import collections, patches

def create_reference_domain_quadrangulation(n_grid_points_x, n_grid_points_y, parametric_space_shape):
  """
  Create a quadrangulation of the reference domain, as needed in the algorithm to create 3D hex meshes. 
  The exact same code is contained in `stl_create_mesh.py` in the function `create_planar_mesh`.
  :param n_grid_points_x: number of grid points in x direction
  :param n_grid_points_y: number of grid points in y direction
  :param parametric_space_shape: which quadrangulation domain and scheme to use: 0 = unit circle, 1 = unit square, 2 = square with adjusted grid points, 3 = unit circle with adjusted grid points
  :return: returns a list of points
  """

  n_grid_points = n_grid_points_x*n_grid_points_y
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
        
        if n_grid_points_x%2 == 1 and i == int(n_grid_points_x/2) and j == int(n_grid_points_y/2):   # center point
          x = 0.
          y = 0.
        
      grid_points_parametric_space[j*n_grid_points_x+i] = np.array([x,y])
      
  return grid_points_parametric_space

if __name__ == "__main__":
  # execute only if run as a script
    
  # define the number of grid points
  n_grid_points_x = 11
  n_grid_points_y = 11

  if len(sys.argv) > 1:
    n_grid_points_x = (int)(sys.argv[1])
    n_grid_points_y = n_grid_points_x
    
  if len(sys.argv) == 1:
    print("usage: plot_quadrangulation_schemes.py [<n_grid_points_x>]")

  print("Number of grid points: {} x {}".format(n_grid_points_x, n_grid_points_y))

  # The following code will plot 2x2 subplots with all four reference domain quadrangulations
  # ------------------------------------------------------------------------------------------

  # define the quadrangulation schemes to use
  # 0 = unit circle, 1 = square, 2 = square with adjusted grid points, 3 = unit circle with adjusted grid points
  parametric_space_shapes = [[1, 2], [0, 3]]

  # define the titles of the four plots
  titles = [["Unit square, scheme 1", "Unit square, scheme 2"], ["Unit circle, scheme 1", "Unit circle, scheme 2"]]

  # create the subplots
  f, ax = plt.subplots(2,2,figsize=(10,10))

  # loop over the subplots
  for axis_j in range(2):
    for axis_i in range(2):
      
      # create the quadrangulation
      parametric_space_shape = parametric_space_shapes[axis_j][axis_i]
      grid_points_parametric_space = create_reference_domain_quadrangulation(n_grid_points_x, n_grid_points_y, parametric_space_shape)
      
      # create auxiliary data structures that are needed for the plots
      patches_parametric = []
      parametric_points = []
      factor = 1.0
      
      # loop over grid points in parametric space
      for (j,y) in enumerate(np.linspace(0.0,1.0,n_grid_points_y)):
        phi = float(y) * (n_grid_points_y-1.0) / n_grid_points_y  * 2.*np.pi
        for (i,x) in enumerate(np.linspace(0.0,1.0,n_grid_points_x)):
      
          [x,y] = grid_points_parametric_space[j*n_grid_points_x+i]
      
          # for debugging create markers at grid points in parametric space

          if parametric_space_shape == 1 or parametric_space_shape == 2 or parametric_space_shape == 3 or parametric_space_shape == 4:  # unit square  
            if i == n_grid_points_x-1 or j == n_grid_points_x-1:
              continue
          if parametric_space_shape == 0:
            if i == n_grid_points_x-1:
              continue
          
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

      # plot the quadrilateral elements
      p = collections.PatchCollection(patches_parametric,edgecolors="k",facecolors="white")
      ax[axis_j,axis_i].add_collection(p)
      ax[axis_j,axis_i].plot([p[0] for p in parametric_points],[p[1] for p in parametric_points], 'ko')
      
      # plot the points as markers, the boundary points are colored
      for j in range(n_grid_points_y):
        for i in range(n_grid_points_x):
          p = grid_points_parametric_space[j*n_grid_points_x+i]
          if j == 0:
            ax[axis_j,axis_i].plot(p[0], p[1], 'ro')
          elif j == n_grid_points_y-1:
            ax[axis_j,axis_i].plot(p[0], p[1], 'go')
          elif i == 0: 
            ax[axis_j,axis_i].plot(p[0], p[1], 'bo')
          elif i == n_grid_points_x-1: 
            ax[axis_j,axis_i].plot(p[0], p[1], 'yo')
          else:
            ax[axis_j,axis_i].plot(p[0], p[1], 'ko')
          
      # set the title of the subplot
      ax[axis_j,axis_i].set_title(titles[axis_j][axis_i])
      
      # set the x and y limits of the subplot
      if parametric_space_shape == 1 or parametric_space_shape == 2:
        ax[axis_j,axis_i].set_xlim(-0.1,1.1)
        ax[axis_j,axis_i].set_ylim(-0.1,1.1)
      else:
        ax[axis_j,axis_i].set_xlim(-1.1,1.1)
        ax[axis_j,axis_i].set_ylim(-1.1,1.1)
      ax[axis_j,axis_i].axis('equal')
    
  # show plot
  plt.show()
  plt.close()
    
