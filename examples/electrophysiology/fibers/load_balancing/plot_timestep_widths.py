#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import matplotlib.pyplot as plt

filename = "build_release/out/log_dt"

if len(sys.argv) > 1:
  filename = sys.argv[1]

print("filename: {}".format(filename))

data = [
  np.genfromtxt(filename+".0.csv", delimiter=';'),
  np.genfromtxt(filename+".1.csv", delimiter=';'),
  np.genfromtxt(filename+".2.csv", delimiter=';'),
  np.genfromtxt(filename+".3.csv", delimiter=';')]
  
n_points = np.size(data[0],0)

print("loaded {} data points".format(n_points))

fig, axs = plt.subplots(4,figsize=(12,9))
for i in range(4):
  axs[i].plot(data[i][:,0], data[i][:,1])
  axs[i].set_yscale('log')
  axs[i].set_title("rank {}".format(i))
  axs[i].set_ylabel('timestep width [ms]')
  
axs[3].set_xlabel('time [ms]')
plt.tight_layout()

plot_filename = filename+".png"
plt.savefig(plot_filename)
print("created plot \"{}\"".format(plot_filename))

plt.show()
