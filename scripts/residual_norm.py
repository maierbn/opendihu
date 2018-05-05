#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script to plot residual norm which is stored in a file.
# Such a file will be created by the nonlinear solver if the "logfile": "filename" option is set.
# The file contains lines "index;number", e.g. 0:2.5698
#
# Arguments: [<residual norm text file>]
#

import sys, os
import matplotlib.pyplot as plt

filenames = ["residual_norm.txt"]
if len(sys.argv) > 1:
  filenames = sys.argv[1:]

print filenames

# prepare plot
fig = plt.figure()

# loop over filenames 
for filename in filenames:

  # read data
  try:
    logfile = open(filename, "r")
  except:
    print("Could not open file \"{}\"".format(filename))
    continue

  xlist = []
  ylist = []

  # loop over lines in file
  parsing_failed = False
  for line in logfile.readlines():
    try:
      entries = line.split(";")
      x_entry = int(entries[0])
      y_entry = float(entries[1])
    except:
      parsing_failed = True
      break
      
    xlist.append(x_entry)
    ylist.append(y_entry)
    
  logfile.close()

  if not parsing_failed:
    plt.plot(xlist, ylist, 'o-', label=os.path.splitext(filename)[0])

plt.gca().set_ylim(bottom=0)

plt.grid(which='major')
plt.legend(loc='best')

# show plot window
plt.show()
