#!/usr/bin/python 

import numpy as np
from matplotlib import pyplot as plt

def linear_lagrange0(xi):
  return 1-xi
  
def linear_lagrange1(xi):
  return xi

def quadratic_lagrange0(xi):
  return 2*xi*xi - 3*xi + 1
  
def quadratic_lagrange1(xi):
  return -4*xi*xi + 4*xi

def quadratic_lagrange2(xi):
  return 2*xi*xi - xi
  
def hermite0(xi):
  return 1 - 3*xi*xi + 2*xi*xi*xi
  
def hermite1(xi):
  return xi * (xi-1) * (xi-1)

def hermite2(xi):
  return xi*xi * (3 - 2*xi)

def hermite3(xi):
  return xi*xi * (xi-1)

xlist = np.linspace(0,1,100)

# as plot
plt.rcParams.update({'font.size': 20, 'lines.linewidth': 8, 'lines.markersize':8})
plt.figure()
plt.plot(xlist, [linear_lagrange0(x) for x in xlist], label=r"$\phi_1$")
plt.plot(xlist, [linear_lagrange1(x) for x in xlist], label=r"$\phi_2$")
plt.legend(loc="best")
plt.xlabel(r"$\xi$")
plt.savefig("lagrange.pdf")

# as geometry
# linear lagrange
plt.figure(figsize=(5,5))
plt.plot(xlist, [linear_lagrange0(x) for x in xlist], label=r"$\phi_1$")
plt.plot(xlist, [linear_lagrange1(x) for x in xlist], label=r"$\phi_2$")
plt.plot(xlist, [0 for x in xlist], color='k', lw=4)
#plt.legend(loc="best")
plt.xlabel(r"$\xi$")
plt.axis('off')
plt.savefig("lagrange1.pdf")

# quadratic lagrange
plt.figure(figsize=(5,5))
plt.plot(xlist, [quadratic_lagrange0(x) for x in xlist], label=r"$\phi_1$")
plt.plot(xlist, [quadratic_lagrange1(x) for x in xlist], label=r"$\phi_2$")
plt.plot(xlist, [quadratic_lagrange2(x) for x in xlist], label=r"$\phi_3$")
plt.plot(xlist, [0 for x in xlist], color='k', lw=4)
#plt.legend(loc="best")
plt.xlabel(r"$\xi$")
plt.axis('off')
plt.savefig("lagrange2.pdf")

# hermite
plt.figure(figsize=(5,5))
plt.plot(xlist, [hermite0(x) for x in xlist], label=r"$\phi_1$")
plt.plot(xlist, [hermite1(x) for x in xlist], label=r"$\phi_1$")
plt.plot(xlist, [hermite2(x) for x in xlist], label=r"$\phi_1$")
plt.plot(xlist, [hermite3(x) for x in xlist], label=r"$\phi_1$")
plt.plot(xlist, [0 for x in xlist], color='k', lw=4)
#plt.legend(loc="best")
plt.xlabel(r"$\xi$")
plt.axis('off')
plt.savefig("hermite.pdf")




plt.show()
