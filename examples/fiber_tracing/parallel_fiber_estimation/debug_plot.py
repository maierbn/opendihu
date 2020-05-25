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

if False:
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

# case 2:
p0=np.array([4.5147446983160675,14.269668406973052,1.6938222646713257]); p1=np.array([4.600984389910788,14.266659718833726,1.6938222646713257]); p2=np.array([4.625363384043174,14.360415541369274,1.693822264671326]); p3=np.array([4.806398144347733,14.420189111046364,1.6938222646713255]) # 1+1+1+1<3

points=[p0,p1,p2,p3]
x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]
# points
ax.scatter(x_list[0],y_list[0],z_list[0],c='r',s=30)
ax.scatter(x_list, y_list, z_list, c='b',s=20)

# 2 3
# 0 1

# lines
line_indices = [[0,1],[1,3],[3,2],[2,0]]
x_lines_list = []
y_lines_list = []
z_lines_list = []
for [i1,i2] in line_indices:
  x_lines_list += [x_list[i1], x_list[i2], np.inf]
  y_lines_list += [y_list[i1], y_list[i2], np.inf]
  z_lines_list += [z_list[i1], z_list[i2], np.inf]
ax.plot(x_lines_list, y_lines_list, z_lines_list)

p0=np.array([4.625363384043174,14.360415541369274,1.693822264671326]); p1=np.array([4.806398144347733,14.420189111046364,1.6938222646713255]); p2=np.array([4.756026721500549,14.457674405235165,1.6938222646713257]); p3=np.array([4.907955520751644,14.499204556275572,1.693822264671326]) # 1+1+1+1<3

points=[p0,p1,p2,p3]
x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]
# points
ax.scatter(x_list[0],y_list[0],z_list[0],c='r',s=30)
ax.scatter(x_list, y_list, z_list, c='b',s=20)

# 2 3
# 0 1

# lines
line_indices = [[0,1],[1,3],[3,2],[2,0]]
x_lines_list = []
y_lines_list = []
z_lines_list = []
for [i1,i2] in line_indices:
  x_lines_list += [x_list[i1], x_list[i2], np.inf]
  y_lines_list += [y_list[i1], y_list[i2], np.inf]
  z_lines_list += [z_list[i1], z_list[i2], np.inf]
ax.plot(x_lines_list, y_lines_list, z_lines_list)

p0=np.array([4.756026721500549,14.457674405235165,1.6938222646713257]); p1=np.array([4.907955520751644,14.499204556275572,1.693822264671326]); p2=np.array([4.8907033768309,14.54407974090634,1.6938222646713257]); p3=np.array([5.041876610629544,14.600674934833364,1.6938222646713257]) # 1+1+1+1<3

points=[p0,p1,p2,p3]
x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]
# points
ax.scatter(x_list[0],y_list[0],z_list[0],c='r',s=30)
ax.scatter(x_list, y_list, z_list, c='b',s=20)

# 2 3
# 0 1

# lines
line_indices = [[0,1],[1,3],[3,2],[2,0]]
x_lines_list = []
y_lines_list = []
z_lines_list = []
for [i1,i2] in line_indices:
  x_lines_list += [x_list[i1], x_list[i2], np.inf]
  y_lines_list += [y_list[i1], y_list[i2], np.inf]
  z_lines_list += [z_list[i1], z_list[i2], np.inf]
ax.plot(x_lines_list, y_lines_list, z_lines_list)

p0=np.array([4.600984389910788,14.266659718833726,1.6938222646713257]); p1=np.array([4.768504746767351,14.32978488379697,1.6938222646713255]); p2=np.array([4.806398144347733,14.420189111046364,1.6938222646713255]); p3=np.array([4.925207664672388,14.454329371644802,1.6938222646713257]) # 1+1+1+1<3

points=[p0,p1,p2,p3]
x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]
# points
ax.scatter(x_list[0],y_list[0],z_list[0],c='r',s=30)
ax.scatter(x_list, y_list, z_list, c='b',s=20)

# 2 3
# 0 1

# lines
line_indices = [[0,1],[1,3],[3,2],[2,0]]
x_lines_list = []
y_lines_list = []
z_lines_list = []
for [i1,i2] in line_indices:
  x_lines_list += [x_list[i1], x_list[i2], np.inf]
  y_lines_list += [y_list[i1], y_list[i2], np.inf]
  z_lines_list += [z_list[i1], z_list[i2], np.inf]
ax.plot(x_lines_list, y_lines_list, z_lines_list)

p0=np.array([4.806398144347733,14.420189111046364,1.6938222646713255]); p1=np.array([4.925207664672388,14.454329371644802,1.6938222646713257]); p2=np.array([4.907955520751644,14.499204556275572,1.693822264671326]); p3=np.array([5.035636461061356,14.531944854561006,1.6938222646713257]) # 1+1+1+1<3

points=[p0,p1,p2,p3]
x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]
# points
ax.scatter(x_list[0],y_list[0],z_list[0],c='r',s=30)
ax.scatter(x_list, y_list, z_list, c='b',s=20)

# 2 3
# 0 1

# lines
line_indices = [[0,1],[1,3],[3,2],[2,0]]
x_lines_list = []
y_lines_list = []
z_lines_list = []
for [i1,i2] in line_indices:
  x_lines_list += [x_list[i1], x_list[i2], np.inf]
  y_lines_list += [y_list[i1], y_list[i2], np.inf]
  z_lines_list += [z_list[i1], z_list[i2], np.inf]
ax.plot(x_lines_list, y_lines_list, z_lines_list)

p0=np.array([4.907955520751644,14.499204556275572,1.693822264671326]); p1=np.array([5.035636461061356,14.531944854561006,1.6938222646713257]); p2=np.array([5.041876610629544,14.600674934833364,1.6938222646713257]); p3=np.array([5.208228258827136,14.630743910763027,1.6938222646713257]) # 1+1+1+1<3

points=[p0,p1,p2,p3]
x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]
# points
ax.scatter(x_list[0],y_list[0],z_list[0],c='r',s=30)
ax.scatter(x_list, y_list, z_list, c='b',s=20)

# 2 3
# 0 1

# lines
line_indices = [[0,1],[1,3],[3,2],[2,0]]
x_lines_list = []
y_lines_list = []
z_lines_list = []
for [i1,i2] in line_indices:
  x_lines_list += [x_list[i1], x_list[i2], np.inf]
  y_lines_list += [y_list[i1], y_list[i2], np.inf]
  z_lines_list += [z_list[i1], z_list[i2], np.inf]
ax.plot(x_lines_list, y_lines_list, z_lines_list)

p0=np.array([4.768504746767351,14.32978488379697,1.6938222646713255]); p1=np.array([4.9504811412053185,14.38590944028343,1.693822264671326]); p2=np.array([4.925207664672388,14.454329371644802,1.6938222646713257]); p3=np.array([5.1237068447662555,14.42881220686678,1.6938222646713255]) # 0+1+1+1<3

points=[p0,p1,p2,p3]
x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]
# points
ax.scatter(x_list[0],y_list[0],z_list[0],c='r',s=30)
ax.scatter(x_list, y_list, z_list, c='b',s=20)

# 2 3
# 0 1

# lines
line_indices = [[0,1],[1,3],[3,2],[2,0]]
x_lines_list = []
y_lines_list = []
z_lines_list = []
for [i1,i2] in line_indices:
  x_lines_list += [x_list[i1], x_list[i2], np.inf]
  y_lines_list += [y_list[i1], y_list[i2], np.inf]
  z_lines_list += [z_list[i1], z_list[i2], np.inf]
ax.plot(x_lines_list, y_lines_list, z_lines_list)

p0=np.array([4.925207664672388,14.454329371644802,1.6938222646713257]); p1=np.array([5.1237068447662555,14.42881220686678,1.6938222646713255]); p2=np.array([5.035636461061356,14.531944854561006,1.6938222646713257]); p3=np.array([5.281527502142019,14.524123761129669,1.6938222646713257]) # 1+1+1+1<3

points=[p0,p1,p2,p3]
x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]
# points
ax.scatter(x_list[0],y_list[0],z_list[0],c='r',s=30)
ax.scatter(x_list, y_list, z_list, c='b',s=20)

# 2 3
# 0 1

# lines
line_indices = [[0,1],[1,3],[3,2],[2,0]]
x_lines_list = []
y_lines_list = []
z_lines_list = []
for [i1,i2] in line_indices:
  x_lines_list += [x_list[i1], x_list[i2], np.inf]
  y_lines_list += [y_list[i1], y_list[i2], np.inf]
  z_lines_list += [z_list[i1], z_list[i2], np.inf]
ax.plot(x_lines_list, y_lines_list, z_lines_list)

p0=np.array([5.035636461061356,14.531944854561006,1.6938222646713257]); p1=np.array([5.281527502142019,14.524123761129669,1.6938222646713257]); p2=np.array([5.208228258827136,14.630743910763027,1.6938222646713257]); p3=np.array([5.363735564936542,14.618705477206964,1.6938222646713255]) # 1+1+1+1<3

points=[p0,p1,p2,p3]
x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]
# points
ax.scatter(x_list[0],y_list[0],z_list[0],c='r',s=30)
ax.scatter(x_list, y_list, z_list, c='b',s=20)

# 2 3
# 0 1

# lines
line_indices = [[0,1],[1,3],[3,2],[2,0]]
x_lines_list = []
y_lines_list = []
z_lines_list = []
for [i1,i2] in line_indices:
  x_lines_list += [x_list[i1], x_list[i2], np.inf]
  y_lines_list += [y_list[i1], y_list[i2], np.inf]
  z_lines_list += [z_list[i1], z_list[i2], np.inf]
ax.plot(x_lines_list, y_lines_list, z_lines_list)

plt.show()
