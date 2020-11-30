#!/usr/bin/python3

# Generate a firing times file that contains the timesteps when motor units are stimulated.
# The actual firing probabilities, i.e. level of overall activation has to be adjusted in the script itself
# by setting the variables min_value and max_value.
#
# usage: ./generate_firing_times <output filename> <number of MUs>

import sys
import random
import time
import numpy as np

if len(sys.argv) < 2:
  print("usage: ./generate_firing_times <output filename> <number of MUs> [<number of timesteps>]")
  sys.exit(0)

# parse command line arguments  
output_filename = sys.argv[1]
n_motor_units = (int)(sys.argv[2])

n_timesteps = 100
if len(sys.argv) >= 4:
  n_timesteps = (int)(sys.argv[3])

# determine firing probabilities for motor units
probabilities = [0.5 for _ in range(n_motor_units)]
for mu_no in range(n_motor_units):
  
  min_value = 0.01      # firing probability for the least active motor unit (last motor unit)
  max_value = 0.4       # firing probability for the most active motor unit (MU no. 1)

  # ansatz value(i) = c1 + c2*exp(i),
  # value(0) = min = c1 + c2  =>  c1 = min - c2
  # value(n-1) = max = min - c2 + c2*exp(n-1)  =>  max = min + c2*(exp(n-1) - 1)  =>  c2 = (max - min) / (exp(n-1) - 1)
  c2 = (max_value - min_value) / (1.02**(n_motor_units-1) - 1)
  c1 = min_value - c2
  probability = c1 + c2*1.02**(n_motor_units-1 - mu_no)
  
  probabilities[mu_no] = probability
  print("  MU {} has firing probability {}".format(mu_no, probability))
  
# write file
with open(output_filename, "w") as f:
  
  for timestep_no in range(n_timesteps):
    
    values = []
    for mu_no in range(n_motor_units):
      
      if random.random() > probabilities[mu_no]:
        values.append("0")
      else:
        values.append("1")
        
    f.write(" ".join(values)+"\n")
    
print("File \"{}\" written.".format(output_filename))
