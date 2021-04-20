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
with open('build_release/result.csv') as csvfile:
  reader = csv.reader(csvfile, delimiter=',')
  
  for scenario_name,strain,stress in reader:
    strain = float(strain)
    stress = float(stress)
    
    # if value is invalid
    if abs(strain) > 4 or abs(stress) > 1000:
      continue
      
    a.append((scenario_name,strain,stress))

# sort entries
a = sorted(a)
data = {}
for scenario_name,strain,stress in a:
  if scenario_name not in data:
    data[scenario_name] = {"strain": [], "stress": []}
  data[scenario_name]["strain"].append(strain)
  data[scenario_name]["stress"].append(stress)

# loop over scenario names
labels = [
  "Fully incompressible (OpenDiHu)",
  "Nearly incompressible (OpenDiHu)",
  "Nearly incompressible, decoupled (OpenDiHu)",
  "Nearly incompressible (FEBio)",
]
scenarios = [
  "incompressible_mooney_rivlin",
  "nearly_incompressible_mooney_rivlin",
  "nearly_incompressible_mooney_rivlin_decoupled",
  "nearly_incompressible_mooney_rivlin_febio",
]
markers = [
  "o-",
  "+-",
  "o-",
  "+-",
]

# define global plotting parameters
plt.rcParams.update({'font.size': 16})
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 10

fig, ax = plt.subplots(figsize=(12,6))

for i,scenario_name in enumerate(scenarios):
  
  ax.plot(data[scenario_name]["strain"], data[scenario_name]["stress"], markers[i], markerfacecolor='none', markeredgewidth=2, label=labels[i])
  
plt.ylabel('stress [$N/cm^2$]')
plt.xlabel('strain [-]')
plt.grid(which='major')
plt.subplots_adjust(right=0.5)
ax.legend(bbox_to_anchor=(1.0, 1.0),frameon=False,borderaxespad=0)
plt.savefig("validation_tensile_test.pdf")
plt.show()
