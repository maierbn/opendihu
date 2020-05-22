#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script to plot result of tensile tests, different numbers of elements
#

import sys
import numpy as np
import matplotlib.pyplot as plt
import csv

# parse all data from result.csv
data = {}
with open('build_release/result_convergence.csv') as csvfile:
  reader = csv.reader(csvfile, delimiter=',')
  
  for scenario_name,nx,strain,stress in reader:
    nx = int(nx)
    strain = float(strain)
    stress = float(stress)
    
    # if value is invalid
    if abs(strain) > 2 or abs(stress) > 10:
      continue
    
    if scenario_name not in data:
      data[scenario_name] = {"nx": [], "strain": [], "stress": []}
    data[scenario_name]["nx"].append(nx)
    data[scenario_name]["strain"].append(strain)
    data[scenario_name]["stress"].append(stress)
    
# loop over scenario names
for scenario_name in data.keys():
  
  plt.plot(data[scenario_name]["nx"], data[scenario_name]["stress"], 'o-', label="{} stress".format(scenario_name))
  plt.plot(data[scenario_name]["nx"], data[scenario_name]["strain"], '+-', label="{} strain".format(scenario_name))
  
plt.ylabel('stress, strain [$N/cm^2$], [-]')
plt.xlabel('nx [-]')
#plt.xscale('log',basex=2)
#plt.yscale('log',basey=2)
plt.grid(which='major')
plt.legend(loc='best')
plt.show()
