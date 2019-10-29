#! /usr/bin/python
# generate an exponentially distributed fiber distribution
# usage: ./generate_fiber_distribution <output filename> <number of MUs> [<n fibers>]

import sys
import random
import numpy as np
import scipy
import scipy.integrate
import scipy.signal
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D

if len(sys.argv) < 3:
  print("usage: ./generate_fiber_distribution <output filename> <number of MUs> [<mode> [<n fibers>]] ")
  print("  mode: 0 = random assignment of fibers to MUs (default)")
  print("        1 = MUs are centered at different locations") 
  sys.exit(0)

# parse command line arguments  
output_filename = sys.argv[1]
n_motor_units = (int)(sys.argv[2])

mode = 0
if len(sys.argv) > 3:
  mode = (int)(sys.argv[3])
  if mode not in [0,1]:
    print("mode was {}, set to 0".format(mode))
    mode = 0
    
n_fibers = 10000
if len(sys.argv) > 4:
  n_fibers = (int)(sys.argv[4])
  
  if mode == 1:
    print("Note, n_fibers was given, but mode=1, therefore no n_fibers is needed.")
  
  
print("output_filename: {}, n_motor_units: {}, mode: {}, n_fibers: {}".format(output_filename, n_motor_units, mode, n_fibers))

factor = 1.20

def pdf(x):
  """ the distribution which the number of fibers per MU should follow """
  return factor**x

def cdf(x):
  """ the cumulative distribution function """
  scaling_factor = scipy.integrate.quad(pdf, 1, n_motor_units)[0]
  
  return scipy.integrate.quad(pdf, 1, x)[0] / scaling_factor

def inverse_cdf(y):
  """ the inverse of the cdf(x) function """
  result = scipy.optimize.minimize_scalar(lambda x: np.linalg.norm(cdf(x) - y))
  return result.x
  
ystart = cdf(0.5)
yend = cdf(n_motor_units+0.5)

def draw_sample():
  return int(np.round(inverse_cdf(random.uniform(ystart,yend))))
  
# initialize random to produce deterministic values
random.seed(0)

# random placement of fibers
if mode == 0:

  # output start and end values of cdf, which should be 0 and 1
  x_list = np.linspace(1,n_motor_units,100)
  print("cdf(1): {}, cdf(max): {}".format(cdf(1), cdf(n_motor_units)))

  # plot pdf
  print("test if inverse is correct: 0=",max([inverse_cdf(cdf(x))-x for x in x_list]))
  plt.plot(x_list, [pdf(x) for x in x_list], label="pdf")
  plt.legend()

  # plot cdf and inverse
  plt.figure(2)
  plt.plot(x_list, [cdf(x) for x in x_list], label="cdf")
  y_list = np.linspace(0,1,100)
  plt.plot(y_list, [inverse_cdf(y) for y in y_list], label="inverse cdf")
  plt.xlabel("motor unit no")
  plt.ylabel("probability for a fiber to be in this MU")
  plt.gca().axis('equal')
  plt.grid(which='both')
  #plt.plot(x_list, [inverse_distribution(x) for x in x_list], label="f^-1")
  plt.legend()
  plt.show()
  
  with open(output_filename,"w") as f:
    values = [draw_sample() for i in range(n_fibers)]
    str_values = map(str, values)
    f.write(" ".join(str_values))
  
# centralized placement of fibers
elif mode == 1:
  
  tolerance = 1e-5
  enable_plots = True
  
  print("The number of fibers in x and y direction has to be known:")
  n_fibers_x = input("Please enter n_fibers_x and press Enter: ")
  n_fibers_y = n_fibers_x

  # generate random positions of the MUs
  mu_positions = []
  for mu_no in range(n_motor_units):
    
    # determine range in which the motor unit should be placed
    # add a margin of 10% at the border where the center of the motor unit should not be
    x_start = (int)(np.round(0.1*(n_fibers_x-1)))
    x_end = (int)(np.round(0.9*(n_fibers_x-1)))
    
    x = random.uniform(x_start, x_end)
    y = random.uniform(x_start, x_end)
    
    mu_positions.append(np.array([x,y]))
  
  print("motor unit positions: ")
  print(mu_positions)
  
  def pdf_distance(i,j,x):
    
    current_position = np.array([i,j])
    
    mu_no = (int)(np.round(x))
    if mu_no < 1:
      mu_no = 1
    if mu_no > n_motor_units:
      mu_no = n_motor_units
    mu_position = mu_positions[mu_no-1]
      
    distance = np.linalg.norm(mu_position - current_position)
      
    sigma = 0.1*n_fibers_x
    probability = scipy.stats.norm.pdf(distance,scale=sigma)
    return probability
    
  def cdf_distance(i,j,x):
    scaling_factor = scipy.integrate.quad(lambda x: pdf_distance(i,j,x), 1, n_motor_units, epsabs=tolerance, epsrel=tolerance)[0]
    
    return scipy.integrate.quad(lambda x: pdf_distance(i,j,x), 1, x, epsabs=tolerance, epsrel=tolerance)[0] / scaling_factor
    
  def inverse_cdf_distance(i,j,y):
    
    """ the inverse of the cdf_distance(x) function """
    result = scipy.optimize.minimize_scalar(lambda x: np.linalg.norm(cdf_distance(i,j,x) - y), tol=tolerance)
    return result.x
    
  def pdf_total(i,j,x):
    #return pdf(x)
    #return pdf_distance(i,j,x)
    return pdf(x) * pdf_distance(i,j,x)
    
  def cdf_total(i,j,x):
    scaling_factor = scipy.integrate.quad(lambda x: pdf_total(i,j,x), 1, n_motor_units, epsabs=tolerance, epsrel=tolerance)[0]
    
    return scipy.integrate.quad(lambda x: pdf_total(i,j,x), 1, x, epsabs=tolerance, epsrel=tolerance)[0] / scaling_factor
    
  def inverse_cdf_total(i,j,y):
    
    """ the inverse of the cdf_total function """
    result = scipy.optimize.minimize_scalar(lambda x: np.linalg.norm(cdf_total(i,j,x) - y), tol=tolerance)
    return result.x

  def draw_sample(i,j):
    ystart = cdf_total(i,j,0.5)
    yend = cdf_total(i,j,n_motor_units+0.5)
    
    return int(np.round(inverse_cdf_total(i,j,random.uniform(ystart,yend))))

  # output start and end values of cdf, which should be 0 and 1
  i = (int)(n_fibers_x/2)
  j = (int)(n_fibers_y/2)
  
  n_points_to_plot = 10
  x_list = np.linspace(1,n_motor_units,n_points_to_plot)
  print("at center, (i,j) = ({},{}): cdf_distance(1): {}, cdf_distance(max): {}".format(i,j,cdf_distance(i,j,1), cdf_distance(i,j,n_motor_units)))

  if enable_plots:
    
    # plot pdf
    #print("test if inverse is correct: 0=",max([inverse_cdf_distance(i,j,cdf_distance(i,j,x))-x for x in x_list]))
    scaling_factor = scipy.integrate.quad(lambda x: pdf_distance(i,j,x), 1, n_motor_units, epsabs=tolerance, epsrel=tolerance)[0]
    plt.plot(x_list, [pdf_distance(i,j,x)/scaling_factor for x in x_list], 'o-', label="pdf")
    plt.gca().set_ylim(0,1)
    plt.xticks(range(n_motor_units+1))
    plt.ylabel("probability for a fiber to be in this MU\n (only according to distance)")
    plt.grid(which='both')
    plt.legend()
    plt.savefig("pdf_distance.pdf")

    # plot pdf total
    #print("test if inverse is correct: 0=",max([inverse_cdf_distance(i,j,cdf_distance(i,j,x))-x for x in x_list]))
    scaling_factor = scipy.integrate.quad(lambda x: pdf_total(i,j,x), 1, n_motor_units, epsabs=tolerance, epsrel=tolerance)[0]
    plt.plot(x_list, [pdf_total(i,j,x)/scaling_factor for x in x_list], 'o-', label="pdf")
    plt.gca().set_ylim(0,1)
    plt.xticks(range(n_motor_units+1))
    plt.ylabel("probability for a fiber to be in this MU")
    plt.grid(which='both')
    plt.legend()
    plt.savefig("pdf_total.pdf")

    # plot cdf and inverse
    plt.figure(2)
    plt.plot(x_list, [cdf_distance(i,j,x) for x in x_list], 'o-', label="cdf")
    y_list = np.linspace(0,1,n_points_to_plot)
    plt.plot(y_list, [inverse_cdf_distance(i,j,y) for y in y_list], 'o-', label="inverse cdf")
    plt.xlabel("motor unit no")
    plt.xticks(range(n_motor_units+1))
    plt.gca().axis('equal')
    plt.grid(which='both')
    #plt.plot(x_list, [inverse_distribution(x) for x in x_list], label="f^-1")
    plt.legend()
    plt.savefig("cdf.pdf")
    #plt.show()
  
  # create actual motor unit numbers
  values = []
  # loop over fibers
  for j in range(n_fibers_y):
    for i in range(n_fibers_x):
      mu_no = draw_sample(i,j)
      print("({},{}): {}".format(i,j,mu_no))
      values.append(mu_no)
      
  # write them to file
  with open(output_filename,"w") as f:
    
    str_values = map(str, values)
    f.write(" ".join(str_values))
  
  # plot fibers in 2D plot
  print("n_motor_units: {}, min: {}".format(n_motor_units, (int)(min(values))))

  colors = cm.rainbow(np.linspace(0, 1, n_motor_units))

  index = 0
  point_colors = []
  for j in range(n_fibers_y):
    for i in range(n_fibers_x):
      mu_no = (int)(values[index])
      point_colors.append(colors[mu_no-1,:])
      #plt.plot(i,j,'o',color=colors[mu_no-1,:])
      
      index += 1

  # plot actual center points of MUs
  for mu_no in range(n_motor_units):
    mu_position = mu_positions[mu_no]
    print(mu_position)
    x = mu_position[0]
    y = mu_position[1]
    
    color = colors[mu_no,:]
    print(mu_position,x,y,color)
    plt.plot(x,y, '+', markersize=24,color=color)

  X,Y = np.meshgrid(range(n_fibers_x),range(n_fibers_y))
  plt.scatter(X,Y,color=point_colors,marker="s")

  plt.gca().set_xlim(0,1.1*n_fibers_x)
  plt.gca().set_ylim(0,n_fibers_y+1)

  legend_elements = []

  for mu_no in range(n_motor_units):
    legend_elements.append(Line2D([0], [0], marker="s", color=colors[mu_no], lw=0, label='MU {}'.format(mu_no+1)))
                           
  ax = plt.gca()
  ax.legend(handles=legend_elements, loc='best')

  plt.axis('equal')
  plt.savefig("2d_fiber_distribution_"+output_filename+".pdf")

  plt.show()

