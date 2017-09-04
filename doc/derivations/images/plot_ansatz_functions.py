#!/usr/bin/python 

import numpy as np
from matplotlib import pyplot as plt

def lagrange1(x):
  return 1-x
  
def lagrange2(x):
  return x

xlist = np.linspace(0,1,100)

plt.rcParams.update({'font.size': 20, 'lines.linewidth': 3, 'lines.markersize':8})
plt.figure()
plt.plot(xlist, [lagrange1(x) for x in xlist], label=r"$\phi_1$")
plt.plot(xlist, [lagrange2(x) for x in xlist], label=r"$\phi_2$")
plt.legend(loc="best")
plt.xlabel(r"$\xi$")
plt.savefig("lagrange.pdf")
plt.show()
