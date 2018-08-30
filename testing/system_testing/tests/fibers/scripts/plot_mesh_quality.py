#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt("mesh_quality2.csv", delimiter=";")
print data


plt.rcParams.update({'font.size': 20})
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 8

fig, ax1 = plt.subplots()

ax1.bar(0,data[0][6],color="b")
ax1.set_ylabel('standard deviation [-]', color="b")

ax2 = ax1.twinx()
ax2.bar(1,data[0][7],color="g")
ax2.set_ylabel('duration [s]', color='g')


##
ax1.bar(6,data[1][6],color="b")
ax2.bar(7,data[1][7],color="g")


##
ax1.bar(3,data[2][6],color="b")
ax2.bar(4,data[2][7],color="g")

##
ax1.bar(9,data[3][6],color="b")
ax2.bar(10,data[3][7],color="g")


##
ax1.bar(12,data[4][6],color="b")
ax2.bar(13,data[4][7],color="g")

ax2.set_ylim(0,100);

locs, labels = plt.xticks()           # Get locations and labels

print locs, labels

plt.xticks((1,4,7,10,13), ("Delaunay,\n unit square","Pie chart,\n unit square","Delaunay,\n unit circle","Pie chart,\n unit circle","optimized \n points location"))  # Set locations and labels



plt.grid(which="both")
plt.savefig("result.png")
plt.show()
