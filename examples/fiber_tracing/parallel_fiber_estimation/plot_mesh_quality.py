#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
now = datetime.datetime.now()
print(" ======= plot_mesh_quality.py =======") 
print(now.strftime("%d/%m/%Y %H:%M:%S"))

import sys, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

# load data
if len(sys.argv) == 2:
  filename = sys.argv[1]
else:
  print("usage: <filename>")
  quit()
data = np.genfromtxt(filename, delimiter=";")

# columns:
# 0 scenario name; filename; n valid fibers; variance element lengths; variance angles
# 1 filename; 
# 2 n valid fibers; 
# 3 variance element lengths; 
# 4 variance angles; 

labels = [
  u'0', u'1', u'2'
]
plt.rcParams.update({'font.size': 18})

x = np.arange(len(labels))  # the label locations
width = 0.4  # the width of the bars

color = (0.8,0.2,0.2)
fig, ax = plt.subplots(figsize=(4.5,4))
rects1 = ax.bar(x, list(data[:,4]), width, color=color)
ax.set_ylabel('Variance of angles [rad${}^2$]')
ax.set_xlabel(r'Max. recursion level $\ell_\mathrm{max}$')
ax.set_xticks(x)
ax.set_xticklabels(labels)
plt.grid(which='major', axis='y')
#ax.set_ylim(0,0.8)

#plt.legend([rects1[0]], ['Standard deviation of \nrelative element lengths', 'Duration [s]'], fontsize=16, loc='upper center')
# Add some text for labels, title and custom x-axis tick labels, etc.

fig.tight_layout()
plt.savefig("mesh_quality.pdf")
plt.show()
