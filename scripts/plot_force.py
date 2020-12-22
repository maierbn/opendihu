#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import sys, os
import matplotlib.pyplot as plt
import numpy as np

filename = "muscle_force.csv"

if len(sys.argv) > 1:
  filename = sys.argv[1]

print("filename: {}".format(filename))

# read data
try:
  logfile = open(filename, "r")
except:
  print("Could not open file \"{}\"".format(filename))
  sys.exit(0)

# loop over lines in file
contents = ""
for line in logfile.readlines():
  contents += line
  
pos = contents.rfind("currentTime")
contents = contents[pos:]

with open("a","w") as f:
  f.write(contents)

data = np.genfromtxt("a",delimiter=";", skip_header=1)
print(data)

# 0 = currentTime;
# 1 = forceBottomX;
# 2 = forceBottomY;
# 3 = forceBottomZ;
# 4 = momentBottomX;
# 5 = momentBottomY;
# 6 = momentBottomZ;
# 7 = forceTopX;
# 8 = forceTopY;
# 9 = forceTopZ;
# 10 = momentTopX;
# 11 = momentTopY;
# 12 = momentTopZ

force_bottom = data[:,3]
force_top = data[:,9]

N = 100

force_bottom_moving_avg = np.convolve(force_bottom, np.ones(N)/N, mode='valid')
force_top_moving_avg = np.convolve(force_top, np.ones(N)/N, mode='valid')

fig = plt.figure()
plt.plot(data[N-1:,0], force_top_moving_avg, label="force proximal")
plt.plot(data[N-1:,0], force_bottom_moving_avg, label="force distal")
plt.title("Muscle forces, moving average with $N={}$".format(N))
plt.ylabel("force [N]")
plt.xlabel("time [ms]")
plt.legend()
plt.show()


