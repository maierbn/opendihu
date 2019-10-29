#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D

filename = "MU_fibre_distribution_3780.txt"

print("usage: ./plot_fibre_distribution_2d.py [<filename> [<n_fibers_x> [<n_fibers_y>]]]")

if len(sys.argv) > 1:
  filename = sys.argv[1]

n_fibers_x = 10
if len(sys.argv) > 2:
  n_fibers_x = (int)(sys.argv[2])
  
n_fibers_y = n_fibers_x
if len(sys.argv) > 3:
  n_fibers_y = (int)(sys.argv[2])

print("filename: {}".format(filename))
print("n fibers: {} x {}".format(n_fibers_x, n_fibers_y))

data = np.genfromtxt(filename)

if len(data) != n_fibers_x*n_fibers_y:
  print("Error! File {} contains {} entries, but given number of fibers is {} x {} = {}.".format(filename, len(data), n_fibers_x, n_fibers_y, n_fibers_x*n_fibers_y))
  sys.exit(0)

n_motor_units = (int)(max(data))
print("n_motor_units: {}, min: {}".format(n_motor_units, (int)(min(data))))

colors = cm.rainbow(np.linspace(0, 1, n_motor_units))
print(colors)

index = 0
point_colors = []
for j in range(n_fibers_y):
  for i in range(n_fibers_x):
    mu_no = (int)(data[index])
    point_colors.append(colors[mu_no-1,:])
    #plt.plot(i,j,'o',color=colors[mu_no-1,:])
    
    index += 1

X,Y = np.meshgrid(range(n_fibers_x),range(n_fibers_y))
plt.scatter(X,Y,color=point_colors,marker="s")

plt.gca().set_xlim(0,1.1*n_fibers_x)
plt.gca().set_ylim(0,n_fibers_y+1)


legend_elements = []

for mu_no in range(n_motor_units):
  legend_elements.append(Line2D([0], [0], marker="s", color=colors[mu_no], lw=0, label='MU {}'.format(mu_no+1)))
                         
ax = plt.gca()
ax.legend(handles=legend_elements, loc='best')

plt.axis('equal')
plt.savefig("2d_fiber_distribution_"+filename+".pdf")

plt.show()
