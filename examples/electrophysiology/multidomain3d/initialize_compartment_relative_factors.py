# Multicompartment 3D, biceps
#
import numpy as np
import scipy.stats
import pickle
import sys

def compute_compartment_relative_factors(mesh_node_positions, fiber_data, motor_units):
    
  # list of fibers, fiber = list of points, point = list with 3 coordinate entries
  n_compartments = len(motor_units)

  # create relative factors for compartments
  #if rank_no == 0:
  #  print("determine relative factors for {} motor units:\n{}".format(n_compartments, motor_units))

  # create data structure with 0
  relative_factors = np.zeros((n_compartments, len(mesh_node_positions)))   # each row is one compartment

  # loop over nodes of mesh
  for node_no,node_position in enumerate(mesh_node_positions):
    node_position = np.array(node_position)
    
    # loop over motor units
    for motor_unit_no,motor_unit in enumerate(motor_units):
      
      # find point on fiber that is closest to current node
      fiber_no = motor_unit["fiber_no"]
      if fiber_no >= len(fiber_data):
        print("Error with motor unit {}, only {} fibers available".format(motor_unit, len(fiber_data)))
      else:
        max_distance = None
        for fiber_point in fiber_data[fiber_no]:
          d = np.array(fiber_point) - node_position
          distance = np.inner(d,d)
          if max_distance is None or distance < max_distance:
            max_distance = distance
            #print("node_position {}, fiber_point {}, d={}, |d|={}".format(node_position, fiber_point, d, np.sqrt(distance)))
        
        distance = np.sqrt(max_distance)
        
        
        gaussian = scipy.stats.norm(loc = 0., scale = motor_unit["standard_deviation"])
        value = gaussian.pdf(distance)*motor_unit["maximum"]
        relative_factors[motor_unit_no][node_no] += value
        #print("motor unit {}, fiber {}, distance {}, value {}".format(motor_unit_no, fiber_no, distance, value))
  
  return relative_factors
