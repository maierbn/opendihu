# This is a helper script that sets a lot of the internal variables which are all defined in variables.py

import numpy as np
import scipy
import pickle
import sys,os
import struct
import argparse
import random
import time
sys.path.insert(0, '..')
import variables    # file variables.py


[nx, ny, nz] = [elem + 1 for elem in variables.n_elements_muscle]

#######################################
# position sensor organs in the 3D mesh

# determine (random) positions of muscle spindles in elasticity mesh
_muscle_spindle_node_nos = []
for muscle_spindle_no in range(variables.n_muscle_spindles):
  i = random.randrange(0,nx)
  j = random.randrange(0,ny)
  k = random.randrange(0,nz)

  dof_no_global = k*nx*ny + j*nx + i
  _muscle_spindle_node_nos.append(dof_no_global)
muscle1_spindle_node_nos = _muscle_spindle_node_nos
# the muscle spindle mesh holds muscle spdindels of both muscles
muscle1_spindle_indices = list(range(variables.n_muscle_spindles))




print(f"Muscle 1 indices in spindle    mesh: {muscle1_spindle_indices[0]:3}...{muscle1_spindle_indices[-1]:3}")

