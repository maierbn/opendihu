#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file gets imported by stl_create_mesh.py

import sys, os
import numpy as np
import matplotlib
import copy
import random

havedisplay = False
if not havedisplay:
  #print("use Agg backend")
  matplotlib.use('Agg')
else:
  print("use Tk backend")

import matplotlib.pyplot as plt
from matplotlib import collections, patches

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

def angle_between_points(p,p0,p1):
  # p1
  # p p0
  return np.arctan2(np.linalg.norm(np.cross(p0-p, p1-p)), np.dot(p0-p, p1-p))

def angle_constraint_is_met(p0,p1,p2,p3):
  # p2 p3
  # p0 p1
  
  # variance of side lengths, favours quadrilaterals with same side lengths
  p01 = p1-p0
  p12 = p3-p1
  p23 = p2-p3
  p30 = p0-p2
  
  # angles
  a0 = np.arctan2(np.linalg.norm(np.cross(p01, -p30)), np.dot(p01, -p30))
  a1 = np.arctan2(np.linalg.norm(np.cross(p12, -p01)), np.dot(p12, -p01))
  a2 = np.arctan2(np.linalg.norm(np.cross(p23, -p12)), np.dot(p23, -p12))
  a3 = np.arctan2(np.linalg.norm(np.cross(p30, -p23)), np.dot(p30, -p23))

  angle_constraint = 20/180.*np.pi
  return abs(a0) >= angle_constraint and abs(a1) >= angle_constraint and abs(a2) >= angle_constraint and abs(a3) >= angle_constraint

def angle_constraint_score(p0,p1,p2,p3):
  # p2 p3
  # p0 p1
  
  # variance of side lengths, favours quadrilaterals with same side lengths
  p01 = p1-p0
  p12 = p3-p1
  p23 = p2-p3
  p30 = p0-p2
  
  # angles
  a0 = abs(np.arctan2(np.linalg.norm(np.cross(p01, -p30)), np.dot(p01, -p30)))
  a1 = abs(np.arctan2(np.linalg.norm(np.cross(p12, -p01)), np.dot(p12, -p01)))
  a2 = abs(np.arctan2(np.linalg.norm(np.cross(p23, -p12)), np.dot(p23, -p12)))
  a3 = abs(np.arctan2(np.linalg.norm(np.cross(p30, -p23)), np.dot(p30, -p23)))

  return (min(30/180.*np.pi, a0) + min(30/180.*np.pi, a1) + min(30/180.*np.pi, a2) + min(30/180.*np.pi, a3))

def perform_laplacian_smoothing(grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, point_indices_list, triangle_list, extent_x, extent_y, loop_no, debugging_stl_output, stl_triangle_lists):
  
  if debugging_stl_output:
    [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,\
      markers_grid_points_parametric_space, markers_grid_points_world_space] = stl_triangle_lists
      
  # iterations: ||: from bottom left, from top right, from bottom right, from top left :||
  i_start_list = [1, n_grid_points_x-2, n_grid_points_x-2, 1]
  j_start_list = [1, n_grid_points_y-2, 1, n_grid_points_y-2]
  i_end_list = [n_grid_points_x-1, 0, 0, n_grid_points_x-1]
  j_end_list = [n_grid_points_y-1, 0, n_grid_points_y-1, 0]
  i_increment_list = [1, -1, -1, 1]
  j_increment_list = [1, -1, 1, -1]
  
  # for total number of smoothing steps
  for k in range(5):
    
    # for interior mesh points
    i_start = i_start_list[k%len(i_start_list)]
    j_start = j_start_list[k%len(j_start_list)]
    i_end = i_end_list[k%len(i_end_list)]
    j_end = j_end_list[k%len(j_end_list)]
    i_increment = i_increment_list[k%len(i_increment_list)]
    j_increment = j_increment_list[k%len(j_increment_list)]
    
    for j in range(j_start,j_end,j_increment):
      for i in range(i_start,i_end,i_increment):
    
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
      
        # compute new point position as average of all 4 or 8 neighbors
        #p_changed = 0.125 * (p0+p1+p2+p3+p4+p5+p6+p7)
        p_changed = 2/12*(p1+p3+p5+p7) + 1/12*(p0+p2+p4+p6)
        #p_changed = 0.25 * (p1+p3+p5+p7)
        
        # discard smoothing step if it would lead to an invalid mesh
        if not is_properly_oriented(p0,p1,p7,p_changed) or not is_properly_oriented(p1,p2,p_changed,p3) \
           or not is_properly_oriented(p7,p_changed,p6,p5) or not is_properly_oriented(p_changed,p3,p5,p4):
          continue
          
        # if new quad violates angle constraint, where all angles have to be > 10 degrees
        if (angle_constraint_is_met(p0,p1,p7,p) and not angle_constraint_is_met(p0,p1,p7,p_changed)) \
          or (angle_constraint_is_met(p1,p2,p,p3) and not angle_constraint_is_met(p1,p2,p_changed,p3)) \
          or (angle_constraint_is_met(p7,p,p6,p5) and not angle_constraint_is_met(p7,p_changed,p6,p5)) \
          or (angle_constraint_is_met(p,p3,p5,p4) and not angle_constraint_is_met(p_changed,p3,p5,p4)):
          continue
        
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
        
              if parametric_space_shape == 1 or parametric_space_shape == 2 or parametric_space_shape == 3:  # unit square  
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
         

def resolve_small_angles(grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, point_indices_list, triangle_list, extent_x, extent_y, loop_no):
  
  factor = (extent_x*extent_y)/2 * 5e-3
  any_point_was_changed = False
  
  # repeatedly iterate over all elements until no more fixes could be achieved
  while True:
    changed_a_point = False
    n_unresolved_bad_angles = 0
    n_resolved_bad_angles = 0
    output_fix = False
    
    # loop over all elements
    for i in range(0,n_grid_points_x-1):
      for j in range(0,n_grid_points_y-1):
        # p2 p3
        # p0 p1
    
        p0 = grid_points_world_space_improved[j*n_grid_points_x+i]
        p1 = grid_points_world_space_improved[j*n_grid_points_x+(i+1)%n_grid_points_x]
        p2 = grid_points_world_space_improved[(j+1)%n_grid_points_y*n_grid_points_x+i]
        p3 = grid_points_world_space_improved[(j+1)%n_grid_points_y*n_grid_points_x+(i+1)%n_grid_points_x]

        # debugging output
        p01 = p1-p0
        p12 = p3-p1
        p23 = p2-p3
        p30 = p0-p2
        
        # angles
        a0 = np.arctan2(np.linalg.norm(np.cross(p01, -p30)), np.dot(p01, -p30))
        a1 = np.arctan2(np.linalg.norm(np.cross(p12, -p01)), np.dot(p12, -p01))
        a2 = np.arctan2(np.linalg.norm(np.cross(p23, -p12)), np.dot(p23, -p12))
        a3 = np.arctan2(np.linalg.norm(np.cross(p30, -p23)), np.dot(p30, -p23))
        angle_constraint = 20/180.*np.pi
        
        #if a0 < angle_constraint or a1 < angle_constraint or a2 < angle_constraint or a3 < angle_constraint:
        #  print("({},{}): constraint met={}, angles {:.2f} {:.2f} {:.2f} {:.2f}".format(i,j,angle_constraint_is_met(p0,p1,p2,p3),
        #    a0*180/np.pi,
        #    a1*180/np.pi,
        #    a2*180/np.pi,
        #    a3*180/np.pi))
  
        if not angle_constraint_is_met(p0,p1,p2,p3):
          indices = [(i,j),(i+1,j),(i,j+1),(i+1,j+1)]
          
          # loop over 4 points of element
          for k in range(4):
            (ii,jj) = indices[k]
            
            # do not consider boundary points, they cannot be changed
            if ii <= 0 or jj <= 0 or ii >= n_grid_points_x-1 or jj >= n_grid_points_y-1:
              continue
            
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
              
            n_tries = 0
            size_factor = factor
            p_changed = np.array(p)
            
            initial_score = angle_constraint_score(p0,p1,p7,p_changed) + angle_constraint_score(p1,p2,p_changed,p3) \
              + angle_constraint_score(p7,p_changed,p6,p5) + angle_constraint_score(p_changed,p3,p5,p4)
            
            # while the score is not yet better (and maximum of 50 tries), deflect point p
            while (\
              (angle_constraint_score(p0,p1,p7,p_changed) + angle_constraint_score(p1,p2,p_changed,p3) \
              + angle_constraint_score(p7,p_changed,p6,p5) + angle_constraint_score(p_changed,p3,p5,p4)) <= initial_score \
              and n_tries < 200):
              
              # pseudo-randomly deflect point p
              p_changed = p + np.array([(random.random()-0.5)*size_factor, (random.random()-0.5)*size_factor, 0])
              size_factor *= 1.025    # 1.05**200 = 17292
              
              # exit loop if size factor is too big
              if size_factor > 1e6:
                n_tries = 200
                
              # discard operation step if it would lead to an invalid mesh
              if not is_properly_oriented(p0,p1,p7,p_changed) or not is_properly_oriented(p1,p2,p_changed,p3) \
                 or not is_properly_oriented(p7,p_changed,p6,p5) or not is_properly_oriented(p_changed,p3,p5,p4):
                continue
                
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
                
                current_score=(angle_constraint_score(p0,p1,p7,p_changed) + angle_constraint_score(p1,p2,p_changed,p3) \
                  + angle_constraint_score(p7,p_changed,p6,p5) + angle_constraint_score(p_changed,p3,p5,p4))

                print("({},{}): {}->{}".format(i, j, initial_score*180/np.pi, current_score*180/np.pi))
                        
                # debugging output
                # angles
                a0 = angle_between_points(p_changed,p7,p1)*180/np.pi
                a1 = angle_between_points(p_changed,p1,p3)*180/np.pi
                a2 = angle_between_points(p_changed,p3,p5)*180/np.pi
                a3 = angle_between_points(p_changed,p5,p7)*180/np.pi
                
                plt.savefig("out/loop_no_{}_ij{},{}_out_{}_{}_{:.1f}->{:.1f}_({:.1f},{:.1f},{:.1f},{:.1f}).png".format(loop_no,i,j,k,n_tries,initial_score*180/np.pi,current_score*180/np.pi,a0,a1,a2,a3))
                plt.close()
              n_tries += 1
              
            # discard smoothing step if it would lead to an invalid mesh
            if not is_properly_oriented(p0,p1,p7,p_changed) or not is_properly_oriented(p1,p2,p_changed,p3) \
               or not is_properly_oriented(p7,p_changed,p6,p5) or not is_properly_oriented(p_changed,p3,p5,p4):
              continue
              
            if n_tries < 200:
              if angle_constraint_is_met(p0,p1,p2,p3):
                print("  \033[0;32mSuccessfully resolved bad angle ({},{}) after {} iterations\033[0m".format(i,j,n_tries))
                n_resolved_bad_angles += 1
              else:
                current_score = (angle_constraint_score(p0,p1,p7,p_changed) + angle_constraint_score(p1,p2,p_changed,p3) \
                  + angle_constraint_score(p7,p_changed,p6,p5) + angle_constraint_score(p_changed,p3,p5,p4))
                print("  improvement regarding bad angles ({},{}) after {} iterations, score: {} -> {} (max: 480)".format(i,j,n_tries, initial_score*180/np.pi, current_score*180/np.pi))
                
              # assign newly found point
              grid_points_world_space_improved[jj*n_grid_points_x+ii] = p_changed
              changed_a_point = True
              any_point_was_changed = True
              
            else:
              n_unresolved_bad_angles += 1
            
    # if there was a resolved self_intersection, restart iterations, otherwise we are done
    if not changed_a_point:
      if n_resolved_bad_angles > 0:
        print("  {} elements with angles < 20 degrees have been improved.".format(n_resolved_bad_angles))
      if n_unresolved_bad_angles > 0:
        print("\033[0;31m  {} elements have angles < 20 degrees.    \033[0m".format(n_unresolved_bad_angles))
      break
      
  return any_point_was_changed
       
       
def resolve_self_intersections(grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, extent_x, extent_y, debugging_stl_output):
  
  factor = (extent_x*extent_y)/2 * 5e-3
  
  # try to resolve self-intersecting quadrilaterals
  # --------------------------------------------------
  # repeatedly iterate over all elements until no more fixes could be achieved
  while True:
    changed_a_point = False
    are_all_elements_properly_oriented = True
    n_unresolved_self_intersections = 0
    output_fix = False
    
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
            
            # do not consider boundary points, they cannot be changed
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
                
              # assign newly found point
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
              n_unresolved_self_intersections += 1
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
      #if are_all_elements_properly_oriented:
      #  print("  all elements are properly oriented")
      break
   
    if n_unresolved_self_intersections > 1:
      print("\033[0;31Abort after {} unresolved self-intersections.    \033[0m".format(n_unresolved_self_intersections))
      sys.exit(-1)

  
def fix_and_smooth_mesh(grid_points_world_space, n_grid_points_x, n_grid_points_y, point_indices_list, triangle_list, extent_x, extent_y, loop_no, debugging_stl_output, stl_triangle_lists):
  """
  This is a helper function for stl_create_mesh.create_planar_mesh, defined in create_planar_mesh.py
  Improve the mesh by fixing invalid quadrilaterals, elements with small angles and applying Laplacian smoothing.
  :param grid_points_world_space: the grid points on the muscle slice
  :param n_grid_points_x: number of grid points in x direction of the final quadrangulation
  :param n_grid_points_y: number of grid points in y direction of the final quadrangulation
  
  :param point_indices_list: a list of the indices into the points array for each triangle of the triangulation
  :param triangle_list: the resulting triangles with their points
  :param debugging_stl_output: if list should be filled with STL triangles that can be output to a STL mesh for debugging
  :param stl_triangle_lists: the debugging lists: [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,markers_grid_points_parametric_space, markers_grid_points_world_space]
  """
  debug = False   # enable debugging output
  if debugging_stl_output:
    [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,\
      markers_grid_points_parametric_space, markers_grid_points_world_space] = stl_triangle_lists

  random.seed(0)
  #print("  improving mesh (fixing self-intersections and Laplacian smoothing)")
  
  factor = (extent_x*extent_y)/2 * 5e-3
  grid_points_world_space_improved = copy.deepcopy(grid_points_world_space)
  
  # output grid
  # set the following options to produce more output files for debugging
  output_pre_fix = False
  output_fix = False
  output_post_fix = False
  if output_pre_fix:
    patches_world_improved = []
    
    # loop over grid points in parametric space
    for (jj,y) in enumerate(np.linspace(0.0,1.0,n_grid_points_y)):
      phi = float(y) * (n_grid_points_y-1.0) / n_grid_points_y  * 2.*np.pi
      for (ii,x) in enumerate(np.linspace(0.0,1.0,n_grid_points_x)):
  
        if parametric_space_shape == 1 or parametric_space_shape == 2 or parametric_space_shape == 3:  # unit square  
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
      
  resolve_self_intersections(grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, extent_x, extent_y, debugging_stl_output)
  
  # output grid
  if output_post_fix:
    patches_world_improved = []
    
    # loop over grid points in parametric space
    for (jj,y) in enumerate(np.linspace(0.0,1.0,n_grid_points_y)):
      phi = float(y) * (n_grid_points_y-1.0) / n_grid_points_y  * 2.*np.pi
      for (ii,x) in enumerate(np.linspace(0.0,1.0,n_grid_points_x)):
  
        if parametric_space_shape == 1 or parametric_space_shape == 2 or parametric_space_shape == 3:  # unit square  
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
      
  random.seed(1)
  
  #for i in range(5):
  #  perform_laplacian_smoothing(grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, point_indices_list, triangle_list, extent_x, extent_y, loop_no, debugging_stl_output, stl_triangle_lists)
  #  return grid_points_world_space_improved
  
  for i in range(25):
    # improve point locations by Laplacian smoothing
    # --------------------------------------------------
    perform_laplacian_smoothing(grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, point_indices_list, triangle_list, extent_x, extent_y, loop_no, debugging_stl_output, stl_triangle_lists)
  
    #resolve_self_intersections(grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, extent_x, extent_y, debugging_stl_output)
  
    # try to improve quadrilaterals with too small angles
    # --------------------------------------------------
    any_point_was_changed = resolve_small_angles(grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, point_indices_list, triangle_list, extent_x, extent_y, loop_no)
    
    resolve_self_intersections(grid_points_world_space_improved, n_grid_points_x, n_grid_points_y, extent_x, extent_y, debugging_stl_output)
  
    # if there was no change, do not do smoothing again
    if not any_point_was_changed:
      break
    
  if i > 1:
    print("({} smoothing iterations)".format(i))
  return grid_points_world_space_improved
 
