#! /usr/bin/python
# generate an exponentially distributed fiber distribution
# usage: ./generate_fiber_distribution <output filename> <number of MUs> [<n fibers>]

import sys
import random
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 3:
  print("usage: ./generate_fiber_distribution <output filename> <number of MUs> [<n fibers>]")
  sys.exit(0)
  
output_filename = sys.argv[1]
n_motor_units = (int)(sys.argv[2])

n_fibers = 10000
if len(sys.argv) > 3:
  n_fibers = (int)(sys.argv[3])
  
print("output_filename: {}, n_motor_units: {}, n_fibers: {}".format(output_filename, n_motor_units, n_fibers))

factor = 1.20

def distribution(x):
  return factor**x

def inverse_distribution(x):
  return np.log(x)/np.log(factor)
  
ystart = distribution(0.5)
yend = distribution(n_motor_units+0.5)

x_list = np.linspace(1,n_motor_units,100)
#plt.plot(x_list, [distribution(x) for x in x_list], label="f")
#plt.plot(x_list, [inverse_distribution(x) for x in x_list], label="f^-1")
#plt.legend()
#plt.show()

print "0=",max([distribution(inverse_distribution(x))-x for x in x_list])

with open(output_filename,"w") as f:
  values = [int(np.round(inverse_distribution(random.uniform(ystart,yend)))) for i in range(n_fibers)]
  f.write(" ".join(map(str, values)))
  
  
