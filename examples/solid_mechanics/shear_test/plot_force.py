#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script to plot result of tensile tests, different forces.
#

import sys
import numpy as np
import matplotlib.pyplot as plt
import csv

# parse all data from result.csv
data = {}
with open('build_release/result.csv') as csvfile:
  reader = csv.reader(csvfile, delimiter=',')
  
  for scenario_name,strain,stress in reader:
    strain = float(strain)
    stress = float(stress)
    
    # if value is invalid
    if abs(strain) > 4 or abs(stress) > 1000:
      continue
    
    if scenario_name not in data:
      data[scenario_name] = {"strain": [], "stress": []}
    data[scenario_name]["strain"].append(strain)
    data[scenario_name]["stress"].append(stress)
    
# loop over scenario names
for scenario_name in data.keys():
  
  plt.plot(data[scenario_name]["strain"], data[scenario_name]["stress"], 'o-', label=scenario_name)
  
plt.ylabel('stress [$N/cm^2$]')
plt.xlabel('strain [-]')
plt.grid(which='major')
plt.legend(loc='best')
plt.show()
