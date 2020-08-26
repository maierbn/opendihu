#! /usr/bin/python
# generate an exponentially distributed fiber distribution
# usage: ./generate_fiber_distribution <output filename> <number of MUs> [<n fibers>]

import sys
import random
import time
import numpy as np
import scipy
import scipy.integrate
import scipy.signal
import matplotlib
matplotlib.use('Agg')

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

def pdf_unscaled(x):
  """ the distribution which the number of fibers per MU should follow """
  return factor**x

def pdf(x):
  """ the probability distribution function which the number of fibers per MU should follow """
  scaling_factor = scipy.integrate.quad(pdf_unscaled, 1, n_motor_units)[0]
  return pdf_unscaled(x) / scaling_factor

def cdf(x):
  """ the cumulative distribution function """
  
  return scipy.integrate.quad(pdf, 1, x)[0]

def inverse_cdf(y):
  """ the inverse of the cdf(x) function """
  result = scipy.optimize.minimize_scalar(lambda x: np.linalg.norm(cdf(x) - y))
  return result.x
  
ystart = cdf(0.5)
yend = cdf(n_motor_units+0.5)

def draw_sample():
  return int(np.round(inverse_cdf(random.uniform(ystart,yend))))
  
# initialize random to produce deterministic values
#random.seed(0)

# mode 0: random placement of fibers
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
  
# mode 1: centralized placement of fibers around center of MU
elif mode == 1:
  
  tolerance = 1e-7
  enable_plots = True
  
  print("The number of fibers in x and y direction has to be known:")
  n_fibers_x = input("Please enter n_fibers_x and press Enter: ")
  n_fibers_y = n_fibers_x

  t_start = time.time()

  # generate random positions of the MUs
  mu_positions = []
  
  # determine range in which the motor unit should be placed
  # add a margin of 10% at the border where the center of the motor unit should not be
  x_start = (int)(np.round(0.1*(n_fibers_x-1)))
  x_end = (int)(np.round(0.9*(n_fibers_x-1)))
    
  x = (n_fibers_x-1)/2.
  y = x
  
  # get random positions according to low-discrepancy series
  alpha1 = 0.5545497
  alpha2 = 0.308517
  for mu_no in range(n_motor_units):
    
    x = (x + alpha1) % 1.0
    y = (y + alpha2) % 1.0
    
    #x = random.uniform(x_start, x_end)
    #y = random.uniform(x_start, x_end)
    
    mu_position = np.array([
      x_start + x*(x_end - x_start),
      x_start + y*(x_end - x_start)
    ])
    mu_positions.append(mu_position)

  print("determined MU positions")  

  # plot actual center points of MUs
  if True:
    # plot MU centers in 2D plot
    colors = cm.rainbow(np.linspace(0, 1, n_motor_units))

    plt.figure(0)
    print("motor unit positions: ")
    for mu_no in range(n_motor_units):
      mu_position = mu_positions[mu_no]
      print(mu_position)
      x = mu_position[0]
      y = mu_position[1]
    
      color = colors[mu_no,:]
      plt.plot(x,y, '+', markersize=24,color=color)
  #plt.show()
  
  sigma = 0.1*n_fibers_x
  a_factor = 1/(sigma**2)

  def pdf_distance_unscaled(equalization_factors,i,j,x):
    """ unscaled version of the probability distribution of the motor unit of fiber (i,j), does not sum to 1 """
    
    current_position = np.array([i,j])
    
    mu_no = (int)(np.round(x))
    if mu_no < 1:
      mu_no = 1
    if mu_no > n_motor_units:
      mu_no = n_motor_units
    mu_position = mu_positions[mu_no-1]
      
    # accurate
    if False:
      # determine distance of fiber to center of MU
      distance = np.linalg.norm(mu_position - current_position)
      
      # probability is gaussian radial basis function
      sigma = 0.1*n_fibers_x    # set standard deviation to 10% of the whole domain
      probability = scipy.stats.norm.pdf(distance,scale=sigma)
  
    else:   # fast
      d = mu_position - current_position
      distance = np.inner(d,d)
      probability = min(0,1 - a_factor*distance)    # hat function with width 1/sqrt(a)

    return probability * equalization_factors[mu_no-1]
    
  def pdf_distance(equalization_factors,i,j,x):
    """ the probability distribution of the motor unit of fiber (i,j) """
    
    scaling_factor = sum([pdf_distance_unscaled(equalization_factors,i,j,t) for t in range(1, n_motor_units+1)])
    return pdf_distance_unscaled(equalization_factors,i,j,x) / scaling_factor
    
  def cdf_distance(equalization_factors,i,j,x):
    """ the cdf of pdf_distance """
    
    # discrete case
    x = (int)(np.round(x))
    value = 0
    for t in range(1,x+1):
      value += pdf_distance(equalization_factors,i,j,t)
    return value
    
    # continuous case, not in use here
    scaling_factor = scipy.integrate.quad(lambda x: pdf_distance(i,j,x), 1, n_motor_units, epsabs=tolerance, epsrel=tolerance)[0]
    return scipy.integrate.quad(lambda x: pdf_distance(i,j,x), 1, x, epsabs=tolerance, epsrel=tolerance)[0] / scaling_factor
    
  def inverse_cdf_distance(equalization_factors,i,j,y):
    """ the inverse of the cdf_distance(x) function """
    
    for t in range(1, n_motor_units+1):
      if cdf_distance(equalization_factors,i,j,t) > y:
        break
        
    return t
    
    # continuous case, not in use here
    result = scipy.optimize.minimize_scalar(lambda x: np.linalg.norm(cdf_distance(i,j,x) - y), tol=tolerance)
    return result.x
    
  def draw_sample(i,j):
    """ draw a sample of the random variable, the motor unit for fiber (i,j) """
    ystart = 0
    yend = cdf_distance(equalization_factors,i,j,n_motor_units)
    
    return int(np.round(inverse_cdf_distance(equalization_factors,i,j,random.uniform(ystart,yend))))

  def objective(equalization_factors):
    """ objective function for the optimization of the equalization_factors, how well the total fiber distribution follows the exponential pdf() distribution """
    
    result = 0
    for mu_no in range(1,n_motor_units+1):
      required_propability_per_fiber = pdf(mu_no)
      #required_propability_per_fiber = 1./n_motor_units
      
      expected_value_mu = sum([pdf_distance(equalization_factors,i,j,mu_no) for j in range(n_fibers_y) for i in range(n_fibers_x)])      propability_per_fiber = expected_value_mu / n_fibers_total
    
      result += (required_propability_per_fiber - propability_per_fiber)**2
    
    return result

  # Explanation of equalization_factors:
  # Because the center points of the motor units are random generated, 
  # some MUs will naturally get more fibers, e.g. if they are far away from all other MUs.
  # The equalization_factors change the probability of each MUs to circumvent this problem.
  # They don't equalize the total expected number of fibers for each MU, they even make them follow the wanted exponential distribution.
  
  # initialize factors to account for different 
  equalization_factors = [1 for i in range(n_motor_units)]
  unity_equalization_factors = [1 for i in range(n_motor_units)]
  
  # print expected values of number of fiber per MU, before applying the equalization_factors
  n_fibers_total = n_fibers_x * n_fibers_y
  for mu_no in range(1,n_motor_units+1):
    expected_value_mu = sum([pdf_distance(equalization_factors,i,j,mu_no) for j in range(n_fibers_y) for i in range(n_fibers_x)])
    print("MU {}, initial number of fibers (expected value): {}".format(mu_no, expected_value_mu))
    
  # find equalization_factors that fulfill the required exponential distribution of MU sizes
  print("\nOptimize factors to obtain exponential distribution of MU sizes")
  bounds = [(0,None) for i in range(n_motor_units)]
  result = scipy.optimize.minimize(objective, equalization_factors, bounds=bounds, options={"disp": True, "maxiter": 1e5})
  print("success: {}, status: {}, message: {}, n objective evaluations: {}, n iterations: {}".format(result.success, result.status, result.message, result.nfev, result.nit))
  equalization_factors = result.x
    
  t_end = time.time()
  print("total duration: {}s".format(t_end-t_start))
  
  # print the result
  print("determined equalization_factors:")
  print(equalization_factors)
    
  # print the final expected number of fibers per MU, the real values will be plotted later
  print("\nfinal expected MU sizes:")
  for mu_no in range(1,n_motor_units+1):
    
    expected_value_mu = sum([pdf_distance(equalization_factors,i,j,mu_no) for j in range(n_fibers_y) for i in range(n_fibers_x)])
    print("mu {}, number of fibers (expected value): {}".format(mu_no, expected_value_mu))
    
  # initialize i,j to the center fiber in the muscle
  i = (int)(n_fibers_x/2)
  j = (int)(n_fibers_y/2)
  
  n_points_to_plot = 100
  
  # output start and end values of cdf, which should be 0 and 1
  x_list = np.linspace(1,n_motor_units,n_points_to_plot)
  print("at center, (i,j) = ({},{}): cdf_distance(0): {} (should be 0), cdf_distance(max): {} (should be 1)".
    format(i,j,cdf_distance(equalization_factors,i,j,0), cdf_distance(equalization_factors,i,j,n_motor_units)))

  if enable_plots:
    
    # plot pdf
    #print("test if inverse is correct: 0=",max([inverse_cdf_distance(i,j,cdf_distance(i,j,x))-x for x in x_list]))
    plt.plot(x_list, [pdf_distance(unity_equalization_factors,i,j,x) for x in x_list], '-', label="pdf")
    plt.gca().set_ylim(0,1)
    plt.xticks(range(n_motor_units+1))
    plt.ylabel("probability for a fiber to be in this MU\n (only according to distance)")
    plt.grid(which='both')
    plt.legend()
    plt.savefig("plots/"+output_filename+"_pdf_distance.pdf")

    # plot pdf total
    plt.figure(2)
    #print("test if inverse is correct: 0=",max([inverse_cdf_distance(i,j,cdf_distance(i,j,x))-x for x in x_list]))
    plt.plot(x_list, [pdf_distance(equalization_factors,i,j,x) for x in x_list], '-', label="pdf")
    plt.gca().set_ylim(0,1)
    plt.xticks(range(n_motor_units+1))
    plt.ylabel("probability for a fiber to be in this MU")
    plt.grid(which='both')
    plt.legend()
    plt.savefig("plots/"+output_filename+"_pdf_total.pdf")

    # plot cdf and inverse
    plt.figure(3)
    plt.plot(x_list, [cdf_distance(equalization_factors,i,j,x) for x in x_list], '-', label="cdf")
    y_list = np.linspace(0,1,n_points_to_plot)
    plt.plot(y_list, [inverse_cdf_distance(equalization_factors,i,j,y) for y in y_list], '-', label="inverse cdf")
    plt.xlabel("motor unit no")
    plt.xticks(range(n_motor_units+1))
    plt.gca().axis('equal')
    plt.grid(which='both')
    #plt.plot(x_list, [inverse_distribution(x) for x in x_list], label="f^-1")
    plt.legend()
    plt.savefig("plots/"+output_filename+"_cdf.pdf")
    #plt.show()
  
  # create actual motor unit numbers by sampling according to the pdf
  values = []
  # loop over fibers
  for j in range(n_fibers_y):
    for i in range(n_fibers_x):
      mu_no = draw_sample(i,j)
      if mu_no > n_motor_units:
        mu_no = n_motor_units
      
      values.append(mu_no)
    print("({},{}): {}".format(i,j,mu_no))
      
  # write all motor units to output file
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
  plt.figure(4)
  for mu_no in range(n_motor_units):
    mu_position = mu_positions[mu_no]
    x = mu_position[0]
    y = mu_position[1]
    
    color = colors[mu_no,:]
    print(mu_position,x,y,color)
    plt.plot(x,y, 'x', markersize=24,markeredgewidth=2,color=color)

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
  plt.savefig("plots/"+output_filename+"_2d_fiber_distribution.pdf")

  # plot total distribution of MUs
  factor = 1.20

  plt.figure(5)
  xlist = np.linspace(1,n_motor_units,5*n_motor_units)
  bins = [x-0.5 for x in range(1,n_motor_units+2)]
  numbers,edges = np.histogram(values,bins)
  print("actual MU sizes: {}".format(numbers))
  plt.hist(values,bins=bins, align='mid', rwidth=0.8)
  a = numbers[-1] / (factor**n_motor_units)
  plt.plot(xlist, [factor**x * a for x in xlist], label='${}^x$'.format(factor))
  #plt.plot(xlist, [np.exp(x)/np.exp(n_motor_units)*bins[-1] for x in xlist])
  plt.xlabel("MU no")
  plt.ylabel("count")
  plt.legend()
  plt.savefig("plots/"+output_filename+"_fiber_distribution.pdf")

  plt.show()

