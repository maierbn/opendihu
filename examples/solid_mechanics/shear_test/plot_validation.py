#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Script to plot result of tensile tests, different forces.
#

import sys
import numpy as np
import matplotlib.pyplot as plt
import csv

# parse all data from result.csv
a = []
with open('build_release/result.csv.backup') as csvfile:
  reader = csv.reader(csvfile, delimiter=',')
  
  for scenario_name,strain,stress0,stress1,stress2,stress3,stress4,stress5 in reader:
    strain = float(strain)
    stress0 = float(stress0)
    stress1 = float(stress1)
    stress2 = float(stress2)
    stress3 = float(stress3)
    stress4 = float(stress4)
    stress5 = float(stress5)
    
    # if value is invalid
    if abs(strain) > 4:
      continue
    
    a.append((scenario_name,strain,stress0,stress1,stress2,stress3,stress4,stress5))
    
# sort entries
a = sorted(a)
data = {}
for scenario_name,strain,stress0,stress1,stress2,stress3,stress4,stress5 in a:
  if scenario_name not in data:
    data[scenario_name] = {"strain": [], "stress": []}
  data[scenario_name]["strain"].append(strain)
  data[scenario_name]["stress"].append([stress0,stress1,stress2,stress3,stress4,stress5])

# loop over scenario names
labels = [
  "Fully incompressible (OpenDiHu)",
  "Nearly incompressible (FEBio)",
]
scenarios = [
  "incompressible_mooney_rivlin",
  "nearly_incompressible_mooney_rivlin_febio"
]
markers = [
  "o-",
  "o-",
  "+-",
  "+-"
]

# define global plotting parameters
plt.rcParams.update({'font.size': 16})
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 10

fig, ax = plt.subplots(figsize=(12,6))
colors = []

labels2 = ["$S_{11}$", "$S_{22}$", "$S_{33}$", "$S_{12}$", "$S_{23}$", "$S_{13}$"]
for component_no in range(6):
  
  scenario_name = "incompressible_mooney_rivlin"
  i = 0
  line, = ax.plot(data[scenario_name]["strain"], [stress[component_no] for stress in data[scenario_name]["stress"]], "-", 
     markerfacecolor='none', markeredgewidth=2, label=labels2[component_no])
  
  scenario_name = "nearly_incompressible_mooney_rivlin_febio"
  i = 1
  ax.plot(data[scenario_name]["strain"], [stress[component_no] for stress in data[scenario_name]["stress"]], ":", 
     color=line.get_color(), markerfacecolor='none', markeredgewidth=2)
  
ax.plot([0,0],[0,0],"k-",label="OpenDiHu")
ax.plot([0,0],[0,0],"k:",label="FEBio")
  
plt.ylabel('stress [$N/cm^2$]')
plt.xlabel('strain [-]')
plt.grid(which='major')
plt.subplots_adjust(right=0.53)
ax.legend(bbox_to_anchor=(1.0, 1.0),frameon=False,borderaxespad=0)
plt.savefig("validation_shear_test.pdf")
plt.show()
