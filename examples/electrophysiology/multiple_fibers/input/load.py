#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle

# fibres, 1D lines with 3D coordinates
fibre_file = "laplace3d_structured_linear"
with open(fibre_file, "rb") as f:
  fibres = pickle.load(f)
  
  # format: 
  # list of fibres
  # each fibre is a list of points like [[67.0862, 144.734, 40.4658], [67.0939, 144.736, 40.5989],  ... ]
  
  #print fibres[0]

# --------------------------------
# 3D mesh
mesh_file = "mesh"
with open(mesh_file, 'rb') as f:
  data = pickle.load(f)

  # data is a python dictionary with the following keys:
  print data.keys()
  
  # node_positions: 3D points, e.g. [[65.088264465332031, 146.19198608398438, 40.0], [64.580612182617188, 146.43875122070312, 40.0], ... ]
  # linear_elements: the 8 node indices per element, e.g. [[0, 1, 11, 12, 121, 122, 132, 133], [1, 2, 12, 13, 122, 123, 133, 134], ... ]
  
  # other keys: 'quadratic_elements, 'n_quadratic_elements_per_coordinate_direction', 'seed_points', 'bottom_nodes', 'top_nodes', 'linear_elements', 'n_linear_elements_per_coordinate_direction'
