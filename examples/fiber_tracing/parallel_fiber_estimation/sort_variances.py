#!/bin/python3

import numpy as np

input_filename="build_release/variances.csv"

# parse all data
lines = []
header = ""
with open(input_filename,"r") as f:
  for line in f:
    if '#' in line:
      header = line
    else:
      lines.append(line)

# sort by variance of angles
sorted_lines = sorted(lines, key=lambda line: (float)(line.split(';')[4]))

# output sorted data
with open("{}_sorted.csv".format(input_filename),"w") as f:
  f.write(header)
  for line in sorted_lines:
    f.write(line)
