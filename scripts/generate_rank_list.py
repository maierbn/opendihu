#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generate a list of cpus to be used on Hazel Hen with the fibers_emg example.
# usage: generate_cpu_list <n_subdomains_x> <n_subdomains_y> <n_subdomains_z> <output_filename> 

import sys

n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1
stride = 1    # stride to skip CPUs
output_filename = ""

if len(sys.argv) > 3:
  try:
    n_subdomains_x = (int)(sys.argv[1])
    n_subdomains_y = (int)(sys.argv[2])
    n_subdomains_z = (int)(sys.argv[3])
  except: 
    pass

if len(sys.argv) > 4:
  try:
    stride = (int)(sys.argv[4])
  except:
    pass

if len(sys.argv) > 5:
  output_filename = sys.argv[5]

if output_filename != "":  
  print("n_subdomains: {} {} {}".format(n_subdomains_x, n_subdomains_y, n_subdomains_z))
  print("output_filename: {}".format(output_filename))

# loop over subdomains = ranks
cpu_no = 0
cpu_list = []
for y in range(n_subdomains_y):
  for x in range(n_subdomains_x):
    for z in range(n_subdomains_z):
      cpu_no = z*n_subdomains_y*n_subdomains_x + y*n_subdomains_x + x
      cpu_list.append(cpu_no)
      cpu_no += 1

for i,no in enumerate(list(cpu_list)):
  cpu_list[no] = i*stride

if output_filename != "":
  with open(output_filename, "w") as outfile:
    for i,cpu_no in enumerate(cpu_list):
      if i != 0:
        outfile.write(",")
      outfile.write("{}".format(cpu_no))
    outfile.write("\n")
else:
  for i,cpu_no in enumerate(cpu_list):
    if i != 0:
      sys.stdout.write(",")
    sys.stdout.write("{}".format(cpu_no))
    sys.stdout.flush()
