#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to visualize EMG output files
# arguments: <filename>
#

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

path = "build_release/out/biceps"

# parse file
if len(sys.argv) > 1:
  filenames = sys.argv[1:]
else:
  print("usage: {} <filename>\n if no filename is given, search in subdirectory {}\n\n".format(sys.argv[0],path))
  
  # get all input data in current directory
  files_at_path = os.listdir(path)

  # sort files by number in file name
  filenames = sorted(files_at_path)

  filenames = [os.path.join(path,f) for f in filenames if ".csv" in f and "found" not in f]
filename = filenames[0]

# parse csv file
data = np.genfromtxt(filename, delimiter=';')
t_list = data[:,1]
  
n_points = (int)(data[0][2])

# determine number of electrodes in x direction
# extract all z positions (along muscle)
position_data = data[0][3:3+3*n_points]
electrode_positions = [np.array(position_data[3*i:3*i+3]) for i in range(n_points)]
z_positions = np.array([position_data[3*i+2] for i in range(n_points)])

# compute differences in z position between neighbour points, this will have a jump (high increment) where a new line starts
differences = z_positions[1:] - z_positions[0:-1]

# determine the indices where the increment is higher than 3 times the mean 
jumps = np.where(differences > np.mean(differences)*3)[0]

# number of points in x/y-direction (across muscle) is position of first jump + 1
n_points_xy = jumps[0] + 1
n_points_z = (int)(n_points / n_points_xy)

# plot electrode positions for debugging
if False:
  from mpl_toolkits.mplot3d import Axes3D
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  xpos = [electrode_positions[i][0] for i in range(n_points)]
  ypos = [electrode_positions[i][1] for i in range(n_points)]
  zpos = [electrode_positions[i][2] for i in range(n_points)]
  colors = [(i/n_points,i/n_points,i/n_points) for i in range(n_points)]
  ax.scatter(xpos,ypos,zpos,c=colors)
  plt.show()

# compute average interelectrode distances
# across muscle (xy)
distances_xy = []
for j in range(n_points_z):
  for i in range(n_points_xy-1):
    point_index = j*n_points_xy + i
    
    p0 = electrode_positions[point_index]
    p1 = electrode_positions[point_index+1]
    
    distance = np.linalg.norm(p0-p1)
    distances_xy.append(distance)

# along muscle (z)
distances_z = []
for j in range(n_points_z-1):
  for i in range(n_points_xy):
    point_index = j*n_points_xy + i
    
    p0 = electrode_positions[point_index]
    p1 = electrode_positions[point_index+n_points_xy]
    
    distance = np.linalg.norm(p0-p1)
    distances_z.append(distance)

# compute inter eletrcode distances in cm
inter_electrode_distance_xy = np.mean(distances_xy)
inter_electrode_distance_z = np.mean(distances_z)

print("Input file \"{}\" contains {} x {} = {} points".format(filename, n_points_xy, n_points_z, n_points))
print("Inter-electrode distance: {:.2f} mm x {:.2f} mm".format(10*inter_electrode_distance_xy, 10*inter_electrode_distance_z))

# create plot
n_plots_x = n_points_xy
n_plots_y = n_points_z
fig = plt.figure(figsize=(5,10))

# labels for entire plot
ax = fig.add_subplot(111)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.grid(False)
plt.xlabel("cross-fiber direction")
plt.ylabel("fiber direction")  
plt.title("sEMG for {} x {} electrodes, t: [{}, {}]".format(n_points_xy, n_points_z, t_list[0], t_list[-1]))

# determine minimum and maximum overall values
all_values = data[:,3+n_points*3:]
minimum_value = np.min(all_values)
maximum_value = np.max(all_values)
print("emg value range [mV]: [{},{}]".format(minimum_value, maximum_value))
print("time range [ms]: [{},{}]".format(t_list[0], t_list[-1]))

# loop over sub plots
for j in range(n_plots_y):
  for i in range(n_plots_x):
    index = j*n_plots_x + i
    
    emg_list = data[:,3+n_points*3+index]
    
    ax = fig.add_subplot(n_plots_y,n_plots_x,index+1)   # nrows, ncols, index, plots start top left and increase to the right
    ax.plot(t_list, emg_list)
    ax.set_ylim(minimum_value,maximum_value)
    #ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])

plt.axis('off')
plt.subplots_adjust(hspace=0, wspace=0)
plt.savefig("emg_plot.pdf")
print("Created file \"emg_plot.pdf\".")
plt.show()

