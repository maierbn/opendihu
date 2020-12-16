#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file gets imported by stl_create_mesh.py

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

#np.seterr(all='raise')

def output_debugging_files(grid_points_parametric_space, grid_points_world_space, grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, parametric_space_shape,
                           u, v, points, n_points_per_face, n_additional_points_on_ring, loop_no, show_plot, extent_x, determine_additional_points_on_ring, grid_points_parametric_space_modified, point_indices_list, triangle_list, stl_triangle_lists):
  """
  This is a helper function for stl_create_mesh.create_planar_mesh, defined in create_planar_mesh.py
  Output various plots for debugging.
  :param grid_points_world_space: the grid points on the muscle slice
  :param n_grid_points_x: number of grid points in x direction of the final quadrangulation
  :param n_grid_points_y: number of grid points in y direction of the final quadrangulation
  :param parametric_space_shape: 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
  
  :param point_indices_list: a list of the indices into the points array for each triangle of the triangulation
  :param triangle_list: the resulting triangles with their points
  :param debugging_stl_output: if list should be filled with STL triangles that can be output to a STL mesh for debugging
  :param stl_triangle_lists: the debugging lists: [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,markers_grid_points_parametric_space, markers_grid_points_world_space]
  """
  
  [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,\
    markers_grid_points_parametric_space, markers_grid_points_world_space] = stl_triangle_lists

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
  scale = 10  

  center_point = np.sum(points,axis=0)/len(points)
  z_value = center_point[2]
  
  x_offset = center_point[0] + extent_x*1.5
  y_offset = center_point[1]
  
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

      if parametric_space_shape == 1 or parametric_space_shape == 2 or parametric_space_shape == 3:  # unit square  
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
      
      grid_point_indices_world_space.append([j*n_grid_points_x+i, j*n_grid_points_x+(i+1)%n_grid_points_x, (j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x, (j+1)%n_grid_points_y*n_grid_points_x+i])
      
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
  #ax[1,1].triplot(xw,yw,grid_point_indices_world_space,color='k')
  for quad in grid_point_indices_world_space:
    ax[1,1].plot([xw[quad[0]], xw[quad[1]]], [yw[quad[0]], yw[quad[1]]], '-k')
    ax[1,1].plot([xw[quad[1]], xw[quad[2]]], [yw[quad[1]], yw[quad[2]]], '-k')
    ax[1,1].plot([xw[quad[2]], xw[quad[3]]], [yw[quad[2]], yw[quad[3]]], '-k')
    ax[1,1].plot([xw[quad[3]], xw[quad[0]]], [yw[quad[3]], yw[quad[0]]], '-k')
  
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
  for quad in grid_point_indices_world_space:
    ax[1,2].plot([xw_improved[quad[0]], xw_improved[quad[1]]], [yw_improved[quad[0]], yw_improved[quad[1]]], '-k')
    ax[1,2].plot([xw_improved[quad[1]], xw_improved[quad[2]]], [yw_improved[quad[1]], yw_improved[quad[2]]], '-k')
    ax[1,2].plot([xw_improved[quad[2]], xw_improved[quad[3]]], [yw_improved[quad[2]], yw_improved[quad[3]]], '-k')
    ax[1,2].plot([xw_improved[quad[3]], xw_improved[quad[0]]], [yw_improved[quad[3]], yw_improved[quad[0]]], '-k')
  
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
    
