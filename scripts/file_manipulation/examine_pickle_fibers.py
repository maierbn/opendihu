#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Script to examine generated streamlines, for pickle files.
# For *.bin files, use examine_bin_fibers.py
#
# usage: ./examine_pickle_fibers.py [<filename.pickle> [<no_plots>]]

import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

#fibre_file = "laplace3d_unstructured_quadratic"
fibre_file = "laplace3d_structured_quadratic"
show_plot = True

if len(sys.argv) >= 2:
  fibre_file = sys.argv[1]
  
if len(sys.argv) == 3:
  show_plot = False

with open(fibre_file, "rb") as f:
  fibres = pickle.load(f)
  
  # format: 
  # list of fibres
  # each fibre is a list of points like [[67.0862, 144.734, 40.4658], [67.0939, 144.736, 40.5989],  ... ]
  
fibre_no = 1
fibre = fibres[fibre_no]

# compute distances between points

print("file {}".format(fibre_file))

fig = plt.figure(figsize=[12., 12.]) #Adjusts the aspect ratio and enlarges the figure (text does not enlarge)

ax = fig.add_subplot(111, projection='3d')

for i,streamline in enumerate(fibres):

  print("fibre {} has {} elements".format(i, len(streamline)))
  distances = []
  directions = [[],[],[]]
  min_distance = None
  pos_min_distance = 0
  
  for j in range(len(streamline)-1):
    v = np.array(streamline[j]) - np.array(streamline[j+1])
    distance = np.linalg.norm(v)
    distances.append(distance)
    
    if distance < min_distance or min_distance is None:
      min_distance = distance
      pos_min_distance = j
      
    
    directions[0].append(v[0])
    directions[1].append(v[1])
    directions[2].append(v[2])
  
  print("max distance: {}".format(max(distances)))
  print("min distance: {} at node {}".format(min_distance, pos_min_distance))
  print(distances[pos_min_distance-2:pos_min_distance+3])
  print(streamline[pos_min_distance-2:pos_min_distance+4])
  
  
  print("mean distance: {}".format(np.mean(distances)))
  print("median distance: {}".format(sorted(distances)[(int)(len(distances)/2)]))
  
  ax.plot([p[0] for p in streamline], [p[1] for p in streamline], [p[2] for p in streamline], 'o-')
  ax.set_title("fibre {}".format(i))
  
if show_plot:
  plt.show()
  
  #print(distances)
  
