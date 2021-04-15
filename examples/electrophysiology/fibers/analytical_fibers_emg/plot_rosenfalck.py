#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to plot the rosenfalck function
#

import sys, os
import numpy as np
import matplotlib.pyplot as plt

def g(z):
  """ Rosenfalck phenomenological model """
  if z >= 0:
    return 96 * z**3 * np.exp(-z) - 90
  else:
    return -90
    
# define global plotting parameters
plt.rcParams.update({'font.size': 14})
plt.rcParams['lines.linewidth'] = 2

plt.figure(1)

xlist = np.linspace(0,20,100)
ylist = [g(x) for x in xlist]


plt.plot(xlist, ylist, linewidth=2)
plt.xlabel("z [mm]")
plt.ylabel("Vm [mV]")

plt.savefig("rosenfalck_function.pdf")
plt.show()
