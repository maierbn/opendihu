#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import csv
import collections
import copy
import scipy.spatial
import os
import pickle
import stl_debug_output

import stl
from stl import mesh

import geomdl
from geomdl import NURBS
from geomdl import BSpline
from geomdl import utilities
from geomdl import exchange
from geomdl import construct
from geomdl import linalg
from geomdl import operations
from geomdl import fitting
from geomdl import exchange

def get_u_by_z_value(curve, value):
  """
  Given a B-spline curve c(u) = [x,y,z], determine the parameter value u for which the given z value is reached.
  :param curve: B-spline curve of geomdl package
  :param value: a z value for which the u parameter will be found
  :return: u value such that c(u)[2] = value
  """
  
  debug = False
  
  if debug:
    print("get_u_by_z_value, value={}".format(value))
  
  u = 0.5
  epsilon = 1e-5
  
  # newton search
  increment = 2*epsilon
  while abs(increment) > epsilon:
    #x = curve.evaluate_single(u)
    x,x_prime = curve.derivatives(u, 1)
    x = x[2] - value
    x_prime = x_prime[2]
    increment = -x / x_prime
    #print("u: {}, x: {}, x': {}, increment: {}".format(u, x, x_prime, increment))
    u += increment

  if debug:
    result = curve.evaluate_single(u)[2]
    print("result: u={}, error={}, x={}".format(u,value-result,result))

  return u
  
def get_u_by_xy_value(curve, xyvalue):
  """
  Given a B-spline curve c(u) = [x,y,z], determine the parameter value u for which the point is closest to the given xyvalue
  :param curve: B-spline curve of geomdl package, should be in a x-y-plane
  :param xyvalue: a list [x,y] of a point in the same plane as the curve.
  :return: u value such that c(u) is the closest point to xyvalue
  """
  
  debug = False
  
  if debug:
    print("get_u_by_z_value, xyvalue={}".format(xyvalue))
  
  u = 0.5
  epsilon = 1e-5
  
  # f(s) = sqrt((x(s) - x0)^2 + (y(s) - y0)^2)
  # f'(s) = 1/f(s)*((x(s) - x0)*x'(s) + (y(s) - y0)*y'(s))
  
  # optimization
  def function(u):
    #print("f({})".format(u))
    u = u[0]
    x = curve.evaluate_single(u)
    return np.linalg.norm(np.array(x[:2]) - np.array(xyvalue))
  
  def jacobian(u):
    #print("jacobian({})".format(u))
    u = u[0]
    x,x_prime = curve.derivatives(u, 1)
    f = np.linalg.norm(np.array(x[:2]) - np.array(xyvalue))
    f_prime = 1/f * ((x[0] - xyvalue[0])*x_prime[0] + (x[1] - xyvalue[1])*x_prime[1])
    return np.array([f_prime])
    
  u0 = [0.2]
  result0 = scipy.optimize.minimize(function, u0, jac=jacobian, bounds=[(0,1)])
  
  u0 = [0.8]
  result1 = scipy.optimize.minimize(function, u0, jac=jacobian, bounds=[(0,1)])
  
  if result0["fun"] < result1["fun"]:
    return result0["x"][0]
  else:
    return result1["x"][0]

def get_arc_length(curve, u):
  """
  Determine the curve length of the curve from parameter 0 to u
  :param curve: B-spline curve of geomdl package
  :param u: the final u value until which the curve length is measured
  :return: arc length of curve
  """
  
  #print("get_arc_length({})".format(u))
  
  if u > 1:
    u = 1
  curve.evaluate(start=0.0, stop=u)
  curve_points = curve.evalpts
  
  l = 0
  last_curve_point = curve_points[0]
  for curve_point in curve_points[1:]:
    l += linalg.point_distance(last_curve_point, curve_point)
    last_curve_point = curve_point
  
  #print("  arc_length({})={}".format(u,l))
  
  return l

def get_u_by_arclength(curve, arclength, u_init=0.5):
  """
  Determine the parameter u of the curve such that the arc length of the curve from 0 to u equals arclength.
  :param curve: B-spline curve of geomdl package
  :param arclength: the target length of the curve
  :param u_init: initial guess for u
  :return: u value such that length of curve from 0 to u equals arclength
  """
  
  debug = False
  
  if debug:
    print("get_u_by_arclength, arclength={}".format(arclength))
    
  u = u_init
  epsilon = 1e-5
  
  # newton search
  increment = 2*epsilon
  while abs(increment) > epsilon:
    x = get_arc_length(curve, u) - arclength
    h = 1e-5
    
    if u > h and u < 1-h:
      x_prime = (get_arc_length(curve, u+h) - get_arc_length(curve, u-h)) / (2*h)
    elif u <= h:
      x_prime = (get_arc_length(curve, u+h) - get_arc_length(curve, u)) / h
    else:
      x_prime = (get_arc_length(curve, u) - get_arc_length(curve, u-h)) / h
      
    if debug:
      print("   u: {}, x: {}, x': {}".format(u, x, x_prime))
    increment = -x / x_prime
    u += increment
    if u > 1:
      u = 1
    elif u < 0:
      u = 0
    
  if debug:
    result = get_arc_length(curve, u)
    print("result: u={}, error={}, x={}".format(u,arclength-result,result))
    
  return u

def create_loop(z_value, spline_surface, v_curve, n_points, loop):
  """
  Sample the spline surface at a given z level, result is a loop as list of points
  :param z_value: z level of all extracted points
  :param spline_surface: the spline surface as Surface object of geomdl
  :param v_curve: a b-spline curve in z direction, created by construct.extract_curves(spline_surface)["v"][0]
  :param n_points: number of points that will be sampled
  :param loop: this is the output, this list will contain points
  """
  
  # get the v value for the plane of the curve in z direction
  v = get_u_by_z_value(v_curve, z_value)
  
  # extract u curve in surface at that v value
  [s0, s1] = operations.split_surface_v(spline_surface, v)
  extracted_curves = construct.extract_curves(s1)
  curve = extracted_curves["u"][0]
  
  curve.delta = 0.01
  curve_length = operations.length_curve(curve)
  #curve_length = get_arc_length(curve, 1.0)
  
  u_init = 1e-5
  for s in np.linspace(0,curve_length,n_points+1)[:-1]:
    u_init = get_u_by_arclength(curve,s,u_init)
    point = curve.evaluate_single(u_init)
    point[2] = z_value
    loop.append(point)
  
def create_border_points(spline_surface, bottom_clip, top_clip, n_loops, n_points):
  
  # get curve in z direction
  v_curve = construct.extract_curves(spline_surface)["v"][0]
  
  loops = []
  for z_value in np.linspace(bottom_clip, top_clip, n_loops):
    loop = []
    create_loop(z_value, spline_surface, v_curve, n_points, loop)
    loops.append(loop)
  
  return loops
  
def create_ring_section(spline_surface, start_point, end_point, z_value, n_points, debugging_points=[]):
  """
  Sample points on a curve in the intersection of the x-y-plane z=z_value und the surface mesh, given by spline_surface, from nearest point to start_point to nearest point to end_point.
  The direction is such that the length of the curve is minimal (there are 2 possible orientations cw/ccw).
  :param spline_surface: a spline surface of the geomdl package
  :param start_point: the line starts at the point on the surface with given z_value, that is the nearest to start_point
  :param end_point: the line ends at the point on the surface with given z_value, that is the nearest to end_point
  :param z_value: the z level of the line on the surface
  :param n_points: number of points on the border
  """
  
  debug = False
  
  # get curve in z direction
  v_curve = construct.extract_curves(spline_surface)["v"][0]
  
  # get the v value for the plane of the curve in z direction
  v = get_u_by_z_value(v_curve, z_value)
  
  if debug:
    print("create_ring_section z_value={}, n_points={}, v: {}".format(z_value, n_points, v))
  
  # extract u curve in surface at that v value
  [s0, s1] = operations.split_surface_v(spline_surface, v)
  extracted_curves = construct.extract_curves(s1)
  curve = extracted_curves["u"][0]
  curve_length = operations.length_curve(curve)
  
  # find u value at start_point and end_point
  u_start = get_u_by_xy_value(curve, start_point[0:2])
  u_end = get_u_by_xy_value(curve, end_point[0:2])
  
  if debug:
    print("u in [{}, {}], curve_length: {}".format(u_start, u_end, curve_length))
  
  # sample curve in between
  s_start = get_arc_length(curve, u_start)
  s_end = get_arc_length(curve, u_end)
  
  length_forwards = s_end - s_start
  if length_forwards < 0:
    length_forwards += curve_length
    
  length_backwards = s_start - s_end
  if length_backwards < 0:
    length_backwards += curve_length
  
  if debug:
    print("length forwards: {}, backwards: {}".format(length_forwards, length_backwards))
  
  # if backwards is shorter, change start and end point
  if length_backwards < length_forwards:
    s_start,s_end = s_end,s_start
    if debug:
      print("go backwards, s in [{},{}]".format(s_start,s_end))
  else:
    if debug:
      print("go forwards, s in [{},{}]".format(s_start,s_end))
    
  # forward direction
  points = []
  u = u_start
  
  # if passing border at u=1, u=0 is involved
  if s_end < s_start:
    if debug:
      print("passing border")
    
    section_length = curve_length - s_start + s_end
    interval_length = section_length / (n_points-1)
    
    if debug:
      print("curve_length: {}, section_length: {}, interval_length: {}".format(curve_length, section_length, interval_length))
    
    if debugging_points is not None:
      point = curve.evaluate_single(0.0)
      point[2] = z_value
      debugging_points.append(point)
      point = curve.evaluate_single(1.0)
      point[2] = z_value
      debugging_points.append(point)
    
    s = s_start
    for i in range(n_points):
      if abs(s - s_end) < 1e-7:
        if debug:
          print("reached s_end={}".format(s_end))
        s = s_end
      
      # collect point
      u = get_u_by_arclength(curve, s, u)
      
      if debug:
        print("collect point {} for s={}, u={}".format(i,s,u))
        
      point = curve.evaluate_single(u)
      point[2] = z_value
      points.append(point)
      
      # proceed s to next point
      s += interval_length
      if s > curve_length:
        s -= curve_length
    
  else:
    for s in np.linspace(s_start,s_end,n_points):
      u = get_u_by_arclength(curve,s,u)
      
      if debug:
        print("collect point for s={}, u={}".format(s,u))
      point = curve.evaluate_single(u)
      point[2] = z_value
      points.append(point)
  
  return points
  
# load stl mesh and evaluate 
if __name__ == "__main__":

  if len(sys.argv) < 3:
    print("usage: {} <input filename> <output filename>".format(sys.argv[0]))
    sys.exit(0)

  input_filename = sys.argv[1]
  output_filename = sys.argv[2]
  bottom_clip = 70
  top_clip = 250
  n_loops = 12
  
  n_points_u = 10          # x-y direction (along rings)
  n_points_v = n_loops    # z direction

  
  print("input_filename: {}".format(input_filename))
  print("output_filename: {}".format(output_filename))
  print("bottom_clip: {}".format(bottom_clip))
  print("top_clip: {}".format(top_clip))
  print("input_filename: {}".format(input_filename))
  
  debug = False
  
  # try to load stored surface
  try:
    f = open("surface.pickle.py","rb")
    surface = pickle.load(f)
  except:
    # create surface

    write_output_mesh = False
    rings = stl_create_rings.create_rings(input_filename, bottom_clip, top_clip, n_loops, write_output_mesh)
    [loops, lengths] = stl_create_mesh.rings_to_border_points(rings, n_points_u-1)

    loops = stl_create_mesh.border_point_loops_to_list(loops)

    #print("loops: ",loops)

    # duplicate first point
    new_loops = []
    for loop in loops:
      #new_loop = loop + loop + loop + [loop[0]]
      new_loop = loop + loop + loop + [loop[0]]
      n_points_u = len(new_loop)
      new_loops.append(new_loop)
    loops = new_loops

    # collect points
    evaluated_points = []
    for i in range(n_points_u):
      for j in range(n_points_v):
        evaluated_points.append(loops[j][i])


    print("n_points_u: ",n_points_u)
    print("n_points_v: ",n_points_v)
    print("{} points".format(len(evaluated_points)))
    #print("evaluated_points: ",evaluated_points)
    surface_full = fitting.approximate_surface(evaluated_points, n_points_u, n_points_v, 3, 2, centripetal=True, ctrlpts_size_u=n_points_u-3, ctrlpts_size_v=n_points_v-3)

    [surface0, surface1] = operations.split_surface_u(surface_full, 0.4)
    [surface, surface2] = operations.split_surface_u(surface1, 0.5555)

    # save surface
    with open("surface.pickle.py","wb") as f:
      pickle.dump(surface,f)
  
  # here we have surface
  extracted_curves = construct.extract_curves(surface)
  
  if debug:
    for v_curve in extracted_curves["u"]:
      print("knotvector: {}".format(v_curve.knotvector))
      print("control points: {}".format(v_curve.ctrlpts))
  
  # get the curve in z direction
  v_curve = extracted_curves["v"][0]
  
  if debug:
    print("v_curve knotvector: {}".format(v_curve.knotvector))
    print("v_curve control points: {}".format(v_curve.ctrlpts))
  

  z_level = (bottom_clip + top_clip)/2

  points = []
  n_points = 50
  #create_loop(z_level, surface, v_curve, n_points, points)

    
  start_point=(80,200,146)
  end_point=(92,128,146)
  z_value=146
  n_points=100

  debugging_points=[]
  points = create_ring_section(surface, start_point, end_point, z_value, n_points, debugging_points)

  points.append(start_point)
  points.append(end_point)

  # save resulting points to file
  filename = "points"
  rank_no = 0
  level = -1
  size = 0.1
  stl_debug_output.output_points(filename, rank_no, level, points, size)
  
  stl_debug_output.output_points("debugging", rank_no, level, debugging_points, 0.05)
  
  #sys.exit(0)
  
  for surf in [surface]:
    surf.sample_size = 100
    surf.evaluate()
  
  #surface.vis = VisMPL.VisSurfWireframe()
  if False:
    surface_full.vis = VisPlotly.VisSurface()
    surface_full.render()
    surface.vis = VisPlotly.VisSurface()
    surface.render()
  #VisMPL.save_figure_as(surface.vis,"output.png")
  #exchange.export_stl(surface0, 'surface0.stl')
  #exchange.export_stl(surface1, 'surface1.stl')
  #exchange.export_stl(surface2, 'surface2.stl')
  print("write stl file \"{}\"".format(output_filename))
  exchange.export_stl(surface, output_filename)
  #exchange.export_stl(surface_full, 'surface_full.stl')
  print("File \"{}\" written.".format(output_filename))
  
