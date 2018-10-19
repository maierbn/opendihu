#!../../../../dependencies/python/install/bin/python3
# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# set global parameters for font sizes
plt.rcParams.update({'font.size': 20})
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 8

directory = "build_release"
if len(sys.argv) > 1:
  directory = sys.argv[1]

# plot dt_0D
filename = os.path.join(directory, "dt_0D.csv")
data = np.genfromtxt(filename, skip_header=1, delimiter=";")
xlist = data[:,0]
ylist = data[:,1]
  
# create plot
plt.figure(figsize=(10,8))
plt.grid(which='major')
plt.plot(xlist,ylist,'b+-')
plt.xlabel('dt_0D')
plt.ylabel('relative error')
plt.xscale('log')
plt.yscale('log')

plt.savefig("dt_0D.pdf")

# plot dt_1D
filename = os.path.join(directory, "dt_1D.csv")
data = np.genfromtxt(filename, skip_header=1, delimiter=";")
xlist = data[:,0]
ylist = data[:,1]
  
# create plot
plt.figure(figsize=(10,8))
plt.grid(which='major')
plt.plot(xlist,ylist,'b+-')
plt.xlabel('dt_1D')
plt.ylabel('relative error')
plt.xscale('log')
plt.yscale('log')

plt.savefig("dt_1D.pdf")

# plot dt_3D
filename = os.path.join(directory, "dt_3D.csv")
data = np.genfromtxt(filename, skip_header=1, delimiter=";")
xlist = data[:,0]
ylist = data[:,1]
  
# create plot
plt.figure(figsize=(10,8))
plt.grid(which='major')
plt.plot(xlist,ylist,'b+-')
plt.xlabel('dt_3D')
plt.ylabel('relative error')
plt.xscale('log')
plt.yscale('log')

plt.savefig("dt_3D.pdf")


# plot dt_3Db
filename = os.path.join(directory, "dt_3Db.csv")
data = np.genfromtxt(filename, skip_header=1, delimiter=";")
xlist = data[:,0]
ylist = data[:,1]
  
# create plot
plt.figure(figsize=(10,8))
plt.grid(which='major')
plt.plot(xlist,ylist,'b+-')
plt.xlabel('dt_3D')
plt.ylabel('relative error')
plt.xscale('log')
plt.yscale('log')

plt.savefig("dt_3Db.pdf")

plt.show()
