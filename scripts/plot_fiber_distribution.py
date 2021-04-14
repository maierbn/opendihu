#! /usr/bin/python3
# Plot the fiber to MU assignment that is given in a file, e.g. MU_fibre_distribution_37x37_10.txt
# usage: ./plot_fiber_distribution <filename>
#

import sys
import random
import time
import numpy as np
import scipy
import scipy.integrate
import scipy.signal
import matplotlib
#matplotlib.use('Agg')

import generate_fiber_distribution

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D

if len(sys.argv) <= 3:
  print("usage: ./plot_fiber_distribution <filename> [<n motor units> [<n_fibers_x> [<basis> [<list of mus to plot>]]]]\nIf the parameters are not specified, they are guessed from the data.")
  
  if len(sys.argv) == 1:
    sys.exit(0)
  
basis = None
n_fibers_x = None
n_motor_units = None
mus_to_plot = None
  
# parse command line arguments  
output_filename = sys.argv[1]

if len(sys.argv) > 2:
  n_motor_units = (int)(sys.argv[2])
  
if len(sys.argv) > 3:
  n_fibers_x = (int)(sys.argv[3])
  n_fibers_y = n_fibers_x
  
if len(sys.argv) > 4:
  basis = (float)(sys.argv[4])
  
if len(sys.argv) > 5:
  mus_to_plot = sys.argv[5]
  if "," in mus_to_plot:
    mus_to_plot = list(map(int,mus_to_plot.split(",")))
  else:
    mus_to_plot = [int(mus_to_plot)]

# generate filename with suffix  
output_filename_with_suffix = output_filename
if ".txt" not in output_filename_with_suffix:
  output_filename_with_suffix = "{}.txt".format(output_filename_with_suffix)
  
# load file
f = open(output_filename_with_suffix,"r")
line = f.readline()
mu_nos_for_fibers = line.split(" ")
mu_nos_for_fibers = list(map(int, mu_nos_for_fibers))
  
if n_motor_units is None:
  n_motor_units = max(mu_nos_for_fibers)
  print("Guess number of motor units = {}".format(n_motor_units))
  
if n_fibers_x is None:
  n_fibers = len(mu_nos_for_fibers)
  n_fibers_x = (int)(np.sqrt(n_fibers))
  n_fibers_y = n_fibers_x
  print("Guess n_fibers_x = {}".format(n_fibers_x))
  
if basis is None:
  # compute histogram
  bins = [x-0.5 for x in range(1,n_motor_units+2)]
  numbers,edges = np.histogram(mu_nos_for_fibers,bins)
  xdata = list(range(1,n_motor_units+1))
  
  # fit curve f to data (xdata,numbers)
  import scipy
  def f(x,b):
    scaling_factor_pdf = sum([b**mu_no for mu_no in range(1,n_motor_units+1)])
    return b**x / scaling_factor_pdf * n_fibers_x**2
    
  def df(x,b):
    scaling_factor_pdf = sum([b**x for x in range(1,n_motor_units+1)])
    return x*b**(x-1) / scaling_factor_pdf * n_fibers_x**2
  b,_ = scipy.optimize.curve_fit(f, xdata, numbers, p0=1.1)
  basis = b[0]
  print("Guess basis from data: {}".format(basis))
  
if mus_to_plot is None:
  mus_to_plot = list(range(1,n_motor_units+1))
  
generate_fiber_distribution.basis = basis
generate_fiber_distribution.n_motor_units = n_motor_units
generate_fiber_distribution.scaling_factor_pdf = sum([generate_fiber_distribution.pdf_unscaled(x) for x in range(1,n_motor_units+1)])

mu_position_offset = 0
mu_position_stride = 1
if output_filename_with_suffix[-2:-1] == "_":
  try:
    mu_position_offset = (int)(output_filename_with_suffix[-1])
    print("Assume offset {} for MU positions, i.e., this is one of the 4 subproblems of generation mode 3.".format(mu_position_offset))
    mu_position_stride = 4
  except:
    pass
mu_positions = generate_fiber_distribution.generate_mu_positions(n_motor_units, n_fibers_x, True, output_filename, mu_position_offset, mu_position_stride)

# set global parameters for font sizes
plt.rcParams.update({'font.size': 14})
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 8

print("\n---------------\n Plot file \"{}\"\n {} motor units\n {} x {} fibers\n basis: {}\n MUs to plot: {}".format(output_filename_with_suffix, n_motor_units, n_fibers_x, n_fibers_y, basis, mus_to_plot))
  
print("\nUse the following command to visualize only a subset of the MUs:")
print("  plot_fiber_distribution.py {} {} {} {} 1,2,{}".format(output_filename, n_motor_units, n_fibers_x, round(basis,3), n_motor_units))

# plot fibers in 2D plot
#print("n_motor_units: {}, min: {}".format(n_motor_units, (int)(min(mu_nos_for_fibers))))

colors = cm.rainbow(np.linspace(0, 1, n_motor_units))
color_invalid_mu = np.reshape(np.array((1,1,1,1)), (1,4)) # add white for not assigned fibers
colors = np.concatenate((colors, color_invalid_mu))

index = 0
point_colors = []
for j in range(n_fibers_y):
  for i in range(n_fibers_x):
    mu_no = (int)(mu_nos_for_fibers[index])
    color = colors[mu_no-1,:]
        
    # make mus that are not part of the list of MUs to plot gray
    if mu_no not in mus_to_plot:
        color = (1,1,0.9)

    point_colors.append(color)
    #plt.plot(i,j,'o',color=colors[mu_no-1,:])
    
    index += 1
    
if n_motor_units <= 15:
  fig = plt.figure(5,figsize=(4,4))
else:
  fig = plt.figure(5,figsize=(8,8))
  
# 109 -> 6
# 37  -> 24   (37*24)/n_fibers_x
markersize = min(100,(37*10)/n_fibers_x)

# plot actual center points of MUs
if mu_positions is not None:
  for mu_no in mus_to_plot:
    mu_position = mu_positions[mu_no-1]
    x = mu_position[0]
    y = mu_position[1]
    
    color = colors[mu_no-1,:]
    #print(mu_position,x,y,color)
    plt.plot(x,y, 'x', markersize=24,markeredgewidth=4,color=color)
    plt.plot(x,y, '.', markersize=10,markeredgewidth=1,color="k")

X,Y = np.meshgrid(range(n_fibers_x),range(n_fibers_y))
m = plt.scatter(X,Y,color=point_colors,marker="s",linewidths=None)
m.set_sizes([markersize**2])
           
ax = plt.gca()

legend_elements = []

for mu_no in mus_to_plot:
  legend_elements.append(Line2D([0], [0], marker="s", color=colors[mu_no-1], lw=0, label='MU {}'.format(mu_no)))
        
# reduce the number of legend entries
if len(legend_elements) > 10:
  stride = (int)(len(legend_elements)/10)
  all_legend_elements = list(legend_elements)
  legend_elements = legend_elements[::stride]
  if legend_elements[-1] != all_legend_elements[-1]:
    legend_elements.append(all_legend_elements[-1])
              
ax.legend(handles=legend_elements, bbox_to_anchor=(1.04, 0.5), loc="center left")

ax.set_xlim(0,1.1*n_fibers_x)
ax.set_ylim(0,n_fibers_y+1)
ax.axis('off')
plt.axis('equal')
plt.tight_layout()
plt.savefig("plots/"+output_filename+"_2d_fiber_distribution_.pdf")
plt.savefig("plots/"+output_filename+"_2d_fiber_distribution_.png")

# plot total distribution of MUs

fig = plt.figure(6,figsize=(6.4, 2.8))
bins = [x-0.5 for x in range(1,n_motor_units+2)]
numbers,edges = np.histogram(mu_nos_for_fibers,bins)
print("Smallest MU: {}, Largest MU: {}".format(min(numbers), max(numbers)))

plt.hist(mu_nos_for_fibers,bins=bins, align='mid', rwidth=0.8)
a = numbers[-1] / (basis**n_motor_units)
xlist = np.linspace(1,n_motor_units,5*n_motor_units)

scaling_factor_pdf = sum([basis**mu_no for mu_no in range(1,n_motor_units+1)])
plt.plot(xlist, [basis**x / scaling_factor_pdf * n_fibers_x**2 for x in xlist], lw=4, label='${}^x$'.format(round(basis,3)))
  
#plt.plot(xlist, [np.exp(x)/np.exp(n_motor_units)*bins[-1] for x in xlist])
ax = fig.gca()
ax.set_xlim(1,n_motor_units)
from matplotlib.ticker import MaxNLocator
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.xlabel("Motor Unit Index")
plt.ylabel("Count")
plt.legend()
plt.savefig("plots/"+output_filename+"_fiber_distribution_.pdf")
plt.savefig("plots/"+output_filename+"_fiber_distribution_.png")
plt.show()

print("\nResult plots were written in plot subdirectory as \"plots/"+output_filename+"_*.pdf\".")
#plt.show()

