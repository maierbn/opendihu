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
with open('result.csv') as csvfile:
  reader = csv.reader(csvfile, delimiter=';')
  
  for scenario_name,strain,stress in reader:
    strain = float(strain)
    stress = float(stress)
    
    # if value is invalid
    if abs(stress) < 2:
      continue
    
    if scenario_name not in data:
      data[scenario_name] = {"strain": [], "stress": []}
    data[scenario_name]["strain"].append(strain)
    data[scenario_name]["stress"].append(stress)
    
# loop over scenario names
for scenario_name in data.keys():
  
  plt.plot([1]+data[scenario_name]["strain"], [0]+data[scenario_name]["stress"], 'o-', label=scenario_name)
  
  # compute mean stiffness 
  
  previous_strain = None
  previous_stress = None
  youngs_modulus_list = []
  for strain,stress in zip(data[scenario_name]["strain"],data[scenario_name]["stress"]):
  
    if previous_strain is not None and strain != previous_strain:
      youngs_modulus = (stress - previous_stress) / (strain - previous_strain)
      youngs_modulus_list.append(youngs_modulus)
  
    previous_strain = strain
    previous_stress = stress
  
  E = np.mean(youngs_modulus_list)
  print("{}: Young's modulus: {} N*cm^-2 = {} GPa".format(scenario_name, E, E*1e-5))
plt.ylabel('stress [$N/cm^2$]')
plt.xlabel('strain [-]')
plt.grid(which='major')
plt.legend(loc='best')
plt.show()
