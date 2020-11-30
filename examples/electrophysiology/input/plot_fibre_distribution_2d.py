#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# usage example:
# ./plot_fibre_distribution_2d.py MU_fibre_distribution_combined_67x67_100_2.txt 34 2
# ./plot_fibre_distribution_2d.py MU_fibre_distribution_combined_67x67_100.txt 67
#
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D

filename = "MU_fibre_distribution_3780.txt"

print("usage: ./plot_fibre_distribution_2d.py [<filename> [<n_fibers_x> [<part_no>]]]")

if len(sys.argv) > 1:
  filename = sys.argv[1]

n_fibers_x = 10
if len(sys.argv) > 2:
  n_fibers_x = (int)(sys.argv[2])
  
n_fibers_y = n_fibers_x

print("filename: {}".format(filename))
print("n fibers: {} x {}".format(n_fibers_x, n_fibers_y))

mu_nos_for_fibers = np.genfromtxt(filename)

if len(mu_nos_for_fibers) != n_fibers_x*n_fibers_y:
  print("Error! File {} contains {} entries, but given number of fibers is {} x {} = {}.".format(filename, len(mu_nos_for_fibers), n_fibers_x, n_fibers_y, n_fibers_x*n_fibers_y))
  sys.exit(0)

n_motor_units = (int)(max(mu_nos_for_fibers))
if "sparse" in filename:
  n_motor_units -= 1
print("n_motor_units: {}, min: {}".format(n_motor_units, (int)(min(mu_nos_for_fibers))))

colors = cm.rainbow(np.linspace(0, 1, n_motor_units))
color_invalid_mu = np.reshape(np.array((1,1,1,1)), (1,4))
colors = np.concatenate((colors, color_invalid_mu)) # add white for not assigned fibers

#selected_mus = [11,41,61,81,91]            # use this line instead of the following to only color selected MUs
selected_mus = list(range(1,n_motor_units+1))
  
# plot actual center points of MUs
# generate random positions of the MUs
mu_positions = []

# determine range in which the motor unit should be placed
# add a margin of 10% at the boundary where the center of the motor unit should not be
x_start = (int)(np.round(0.1*(n_fibers_x-1)))
x_end = (int)(np.round(0.9*(n_fibers_x-1)))
  
x = 0.5
y = 0.5
  
# get random positions according to low-discrepancy series
alpha1 = 0.5545497
alpha2 = 0.308517
for mu_no in range(n_motor_units):
  
  x = (mu_no*alpha1) % 1.0
  y = (mu_no*alpha2) % 1.0
  
  mu_position = np.array([
    x_start + x*(x_end - x_start),
    x_start + y*(x_end - x_start)
  ])
  mu_positions.append(mu_position)
  
if "combined" in filename:
  n_motor_units_instance = int(n_motor_units / 4)
  for iteration_no in range(4):
      
    mu_position_offset = iteration_no*n_motor_units_instance
      
    for mu_no in range(mu_position_offset, mu_position_offset+n_motor_units_instance):
      x = (mu_no*alpha1) % 1.0
      y = (mu_no*alpha2) % 1.0
      
      mu_position = np.array([
        x_start + x*(x_end - x_start),
        x_start + y*(x_end - x_start)
      ])
      
      part_mu_no = mu_no-mu_position_offset+1
      index = 4*(part_mu_no-1)+iteration_no
      mu_positions[index] = mu_position
  
if "combined" in filename and len(sys.argv) == 4:
  iteration_no = (int)(sys.argv[3])
  print("Interpreting as result from a combined method for iteration {}".format(iteration_no))
  
  n_motor_units_instance = n_motor_units
  mu_position_offset = iteration_no*n_motor_units_instance

  for mu_no in range(mu_position_offset, mu_position_offset+n_motor_units_instance):
    x = (mu_no*alpha1) % 1.0
    y = (mu_no*alpha2) % 1.0
    
    mu_position = np.array([
      x_start + x*(x_end - x_start),
      x_start + y*(x_end - x_start)
    ])
    
    part_mu_no = mu_no-mu_position_offset+1
    index = part_mu_no-1
    mu_positions[index] = mu_position
  
if mu_positions is not None and False:
  for mu_no in selected_mus:
    if mu_no-1 < 0 or mu_no-1 >= len(mu_positions):
      break
    mu_position = mu_positions[mu_no-1]
    x = mu_position[0]
    y = mu_position[1]
    
    color = colors[mu_no-1,:]
    #print(mu_position,x,y,color)
    plt.plot(x,y, 'x', markersize=24,markeredgewidth=2,color=color)

index = 0
point_colors = []
for j in range(n_fibers_y):
  for i in range(n_fibers_x):
    mu_no = (int)(mu_nos_for_fibers[index])
    color = colors[mu_no-1,:]
      
    # for dataset with high number of motor units only plot some
    if n_motor_units > 20:
      if mu_no not in selected_mus:
        color = (0.8,0.8,0.8)
        
    # no motor unit is white
    if mu_no == n_motor_units+1:
        color = (1.0,1.0,1.0)

    point_colors.append(color)
    
    index += 1

X,Y = np.meshgrid(range(n_fibers_x),range(n_fibers_y))
m = plt.scatter(X,Y,color=point_colors,marker="s")
m.set_sizes([10])
if n_fibers_x < 20:
  m.set_sizes([400])
elif n_fibers_x < 40:
  m.set_sizes([40])
   
ax = plt.gca()

legend_elements = []

for mu_no in selected_mus:
  if mu_no <= n_motor_units:
    legend_elements.append(Line2D([0], [0], marker="s", color=colors[mu_no-1], lw=0, label='MU {}'.format(mu_no)))
        
# reduce the number of legend entries
if True:
  if len(legend_elements) > 10:
    stride = (int)(len(legend_elements)/10)
    all_legend_elements = list(legend_elements)
    legend_elements = legend_elements[::stride]
    if legend_elements[-1] != all_legend_elements[-1]:
      legend_elements.append(all_legend_elements[-1])
                
ax.legend(handles=legend_elements, bbox_to_anchor=(0.85, 0.5), loc="center left")

ax.set_xlim(0,1.1*n_fibers_x)
ax.set_ylim(0,n_fibers_y+1)
ax.axis('off')
plt.axis('equal')
plt.tight_layout()
output_filename = "plots/"+filename+"_2d_fiber_distribution_.pdf"
plt.savefig(output_filename)
print("Saved file \"{}\".".format(output_filename))

# plot total distribution of MUs

fig = plt.figure(6)
bins = [x-0.5 for x in range(1,n_motor_units+2)]
numbers,edges = np.histogram(mu_nos_for_fibers,bins)

plt.hist(mu_nos_for_fibers,bins=bins, align='mid', rwidth=0.8)
ax = fig.gca()
ax.set_xlim(1,n_motor_units)
from matplotlib.ticker import MaxNLocator
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.xlabel("Motor Unit Index")
plt.ylabel("Count")
output_filename = "plots/"+filename+"_fiber_distribution_.pdf"
plt.savefig(output_filename)
print("Saved file \"{}\".".format(output_filename))


plt.show()
