#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to visualize EMG output files
# arguments: <path>
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
  # get all input data in current directory
  files_at_path = os.listdir(path)

  # sort files by number in file name
  filenames = sorted(files_at_path)

  filenames = [os.path.join(path,f) for f in filenames if ".csv" in f and "found" not in f]
filename = filenames[0]

print("file: {}".format(filename))

# parse csv file
data = np.genfromtxt(filename, delimiter=';')
t_list = data[:,1]
  
n_points = (int)(data[0][2])
print("n_points:",n_points)

# determine number of electrodes in x direction
z_positions = np.array([data[0][3+3*i+2] for i in range(n_points)])
differences = z_positions[1:] - z_positions[0:-1]
jumps = np.where(differences > np.mean(differences)*3)[0]
n_points_z = jumps[0] + 1
n_points_xy = (int)(n_points / n_points_z)

n_plots_x = n_points_xy
n_plots_y = n_points_z

print("{} x {} points detected".format(n_points_xy, n_points_z))

fig = plt.figure(figsize=(5,10))

# labels for entire plot
ax = fig.add_subplot(111)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.grid(False)
plt.xlabel("fiber direction")
plt.ylabel("cross-fiber direction")    

# determine minimum and maximum overall values
all_values = data[:,3+n_points*3:]
minimum_value = np.min(all_values)
maximum_value = np.max(all_values)
print("emg value range: [{},{}]".format(minimum_value, maximum_value))
print("time range: [{},{}]".format(t_list[0], t_list[-1]))

# loop over sub plots
for j in range(n_plots_y):
  for i in range(n_plots_x):
    index = j*n_plots_x + i
    
    emg_list = data[:,3+n_points*3+index]
    
    ax = fig.add_subplot(n_plots_x,n_plots_y,index+1)   # nrows, ncolrs, index, plots start top left and increase to the right
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

