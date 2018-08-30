#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script to visualize stl meshes, creates output fig.pdf
#

import sys
import numpy as np

# import needed packages from matplotlib
import matplotlib as mpl
mpl.use('Agg')

from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt

filename = ""
if len(sys.argv) > 1:
  filename = sys.argv[1]

# Create a new plot
fig = plt.figure()
ax = mplot3d.Axes3D(fig)

# Load the STL files and add the vectors to the plot
mesh_object = mesh.Mesh.from_file(filename)
ax.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh_object.vectors))

# Auto scale to the mesh size
scale = mesh_object.points.flatten(-1)
scale *= 0.9
ax.auto_scale_xyz(scale, scale, scale)

# Show the plot to the screen
plt.savefig("fig.pdf")
