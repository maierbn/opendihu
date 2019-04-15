#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import csv
import collections
import copy
from sets import Set
import os


rhs = np.load("p_rhs.npy")
m = np.load("p_stiffness.npy")
solution = np.load("p_solution.npy")
solution_shaped = np.load("p_solution_shaped.npy")

n = rhs.shape[0]

print "n: ", n
print "rhs shape:", rhs.shape
print "m shape: ", m.shape
print "solution shape: ", solution.shape
print "solution_shaped shape: ", solution_shaped.shape

dim = len(solution_shaped.shape)
print "dim: ",dim
f = np.inner(m,solution)

print "f shape: ", f.shape

res = f - rhs
if False:
  for i in range(n):
    print "{}. solution={}, f={}, rhs={}, squared error: {}".format(i, solution[i], f[i], rhs[i], (f[i]-rhs[i])**2)
print "residual norm: ", np.linalg.norm(res)

if dim == 2:

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  x = np.linspace(0,2.0,solution_shaped.shape[0])
  y = np.linspace(0,4.0,solution_shaped.shape[1])
  X, Y = np.meshgrid(x,y)

  print X.shape, Y.shape, solution_shaped.shape
  Z = np.reshape(solution_shaped, [solution_shaped.shape[1],solution_shaped.shape[0]])

  ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=1,rstride=1,cstride=1)

  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')

else:
    
  fig = plt.figure()
  x = np.linspace(0,2.0,solution.shape[0])
  ax = plt.gca()
  ax.plot(x, solution, "o-")

  ax.set_xlabel('X')
  ax.set_ylabel('Y')

plt.show()
