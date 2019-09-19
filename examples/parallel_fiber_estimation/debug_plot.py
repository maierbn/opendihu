#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script to visualize python output files.
#

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

points = [[77.8793,133.399,107.52],[81.4639,129.914,107.52],[74.2754,136.894,107.52],[79.0424,137.149,107.52],[77.8443,133.357,110.48],[81.5631,129.835,110.48],[74.1507,136.931,110.48],[79.4285,136.877,110.48]]

x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]

# points
ax.scatter(x_list[0],y_list[0],z_list[0],c='r',s=30)
ax.scatter(x_list, y_list, z_list, c='b',s=20)

# 2 3
# 0 1
# 
# 6 7
# 4 5

# lines
line_indices = [[0,1],[1,3],[3,2],[2,0],[4,5],[5,7],[7,6],[6,4]]
x_lines_list = []
y_lines_list = []
z_lines_list = []
for [i1,i2] in line_indices:
  x_lines_list += [x_list[i1], x_list[i2], np.inf]
  y_lines_list += [y_list[i1], y_list[i2], np.inf]
  z_lines_list += [z_list[i1], z_list[i2], np.inf]
ax.plot(x_lines_list, y_lines_list, z_lines_list)

plt.show();
