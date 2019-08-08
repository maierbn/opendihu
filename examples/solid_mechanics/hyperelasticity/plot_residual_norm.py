#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script to visualize residual norm
#

import sys
import numpy as np
  
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib import cm
from matplotlib.patches import Polygon

# parse data
indices = []
values = []
with open("log_residual_norm.txt", "r") as f:
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
ax.set_xscale('log')
ax.set_yscale('log')
plt.grid(which='major')
plt.show()
