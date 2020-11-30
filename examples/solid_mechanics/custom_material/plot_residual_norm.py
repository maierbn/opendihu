#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Script to visualize residual norm, that has been written to log file log_residual_norm.txt
#

import sys
import numpy as np
  
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib import cm
from matplotlib.patches import Polygon

filename = "build_release/log_residual_norm.txt"
if len(sys.argv) > 1:
  filename = sys.argv[1]

# parse data
indices = []
values = []
with open(filename, "r") as f:
  lines = f.readlines()
  for line in lines:
    (i,value) = line.split(";")
    indices.append((int)(i)+1)
    values.append((float)(value))

# prepare plot
fig = plt.figure()
ax = plt.gca()
ax.plot(indices, values, 'o-', lw=2)
ax.set_title('convergence of nonlinear solver')
ax.set_xlabel('iteration no.')
ax.set_ylabel('2-norm residual')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
#plt.yscale('log', basey=2)
plt.grid(which='both')
plt.show()
