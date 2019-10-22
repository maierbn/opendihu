#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import sys
import numpy as np
import matplotlib.pyplot as plt

filename = "MU_fibre_distribution_3780.txt"

if len(sys.argv) > 1:
  filename = sys.argv[1]

print("filename: {}".format(filename))

data = np.genfromtxt(filename)

n_motor_units = (int)(max(data))
print("n_motor_units: {}: {}".format(n_motor_units, (int)(min(data))))

bins = [x-0.5 for x in range(1,n_motor_units+2)]
print("bins: {}".format(bins))

factor = 1.20

xlist = np.linspace(1,n_motor_units,5*n_motor_units)
numbers,edges = np.histogram(data,bins)
print(numbers,edges)
plt.hist(data,bins=bins, align='mid', rwidth=0.8)
a = numbers[-1] / (np.exp(n_motor_units))
plt.plot(xlist, [np.exp(x)*a for x in xlist], label='exp(x)')
a = numbers[-1] / (factor**n_motor_units)
plt.plot(xlist, [factor**x * a for x in xlist], label='1.05**x')
#plt.plot(xlist, [np.exp(x)/np.exp(n_motor_units)*bins[-1] for x in xlist])
plt.xlabel("MU no")
plt.ylabel("count")
plt.legend()

plt.savefig("histogram_"+filename+".pdf")
plt.show()
