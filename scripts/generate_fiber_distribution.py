#! /usr/bin/python3
# Generate file that specify assignments of motor units to fibers.
# These files are typically named MU_fibre_distribution_* and stored in examples/electrophysiology/input.
# There are three modes how to generate these distributions. The motor units are always exponentially distributed 
# according to basis^x (by default 1.2^x) with the lower MUs containing smaller number of fibers. The difference between the modes is how
# MU territories are considered.
#
# usage: ./generate_fiber_distribution <output filename> <number of MUs> [<mode> [...]]
# or run without arguments to get help

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
  print("usage: ./generate_fiber_distribution <output filename> <number of MUs> <mode> [<n_fibers_x> [<basis> [...]]] ")
  print("  mode: 0 = random assignment of fibers to MUs (default)")
  print("        1 = MUs are centered at different locations (use this if unsure)")
  print("        2 = like 1, but no scaling on the probabilities is performed,")
  print("            this means numerous fibers will get no MU assigned")
  print("")
  print("  for mode 0:")
  print("       ./generate_fiber_distribution <output filename> <number of MUs> 0 [<n_fibers_x> [<basis>]]\n")
  print("  for mode 1:")
  print("       ./generate_fiber_distribution <output filename> <number of MUs> 1 [<n_fibers_x> [<basis> [<n_max_iterations> [<chunk_size>]]]] ")
  print("       if <n_max_iterations> is set to -1, use pysgpp instead of scipy.optimize (if available)\n")
  print("  for mode 2:")
  print("       ./generate_fiber_distribution <output filename> <number of MUs> 2 [<n_fibers_x> [<basis>]] \n")
  
  sys.exit(0)

# parse command line arguments  
output_filename = sys.argv[1]
n_motor_units = (int)(sys.argv[2])

# parse mode
mode = 0
if len(sys.argv) > 3:
  mode = (int)(sys.argv[3])
  if mode not in [0,1,2]:
    print("mode was {}, set to 0".format(mode))
    mode = 0
    
basis = 1.20

def pdf_unscaled(mu_no):
  """ the distribution which the number of fibers per MU should follow """
  return basis**mu_no

# the pdf scaling factor does only depend on n_motor_units, compute it globally
#scaling_factor_pdf = scipy.integrate.quad(pdf_unscaled, 0.5, n_motor_units+0.5)[0]

scaling_factor_pdf = sum([pdf_unscaled(x) for x in range(1,n_motor_units+1)])

def pdf(mu_no):
  """ the probability distribution function which the number of fibers per MU should follow """
  return pdf_unscaled(mu_no) / scaling_factor_pdf

def cdf(mu_no):
  """ the cumulative distribution function """
  
  return scipy.integrate.quad(pdf, 1, mu_no)[0]

def inverse_cdf(y):
  """ the inverse of the cdf(x) function """
  result = scipy.optimize.minimize_scalar(lambda x: np.linalg.norm(cdf(x) - y))
  return result.x
  
ystart = cdf(0.5)
yend = cdf(n_motor_units+0.5)

def draw_sample():
  return int(np.round(inverse_cdf(random.uniform(ystart,yend))))

def pdf_distance_unscaled(equalization_factors,i,j,mu_no):
  """ unscaled version of the probability distribution of the motor unit of fiber (i,j), does not sum to 1 """
  
  current_position = np.array([i,j])
  
  if mu_no < 1:
    print("error, mu_no is {}".format(mu_no))
    mu_no = 1
  if mu_no > n_motor_units:
    print("error, mu_no is {} > n_motor_units={}".format(mu_no,n_motor_units))
    mu_no = n_motor_units
  mu_position = mu_positions[mu_no-1]
    
  # accurate, use Gaussian RBF
  if use_accurate_radial_basis_function:
    # determine distance of fiber to center of MU
    distance = np.linalg.norm(mu_position - current_position)
    
    # probability is gaussian radial basis function
    sigma = 0.1*n_fibers_x    # set standard deviation to 10% of the whole domain
    probability = scipy.stats.norm.pdf(distance,scale=sigma)

  else:   # fast, use quadratic function
    d = mu_position - current_position
    distance = np.inner(d,d)
    probability = 1/(1 + a_factor*distance)    # quadratic function with width 1/sqrt(a)

  return probability * equalization_factors[mu_no-1]
  
# initialize random to produce deterministic values
#random.seed(0)

# mode 0: random placement of fibers
if mode == 0:

  # parse n_fibers, total number of fibers
  
  n_fibers_x = -1
  if len(sys.argv) > 4:
    n_fibers_x = (int)(sys.argv[4])
  else:
    print("The number of fibers in x and y direction has to be known:")
    n_fibers_x = (int)(input("Please enter n_fibers_x and press Enter: "))
  
  if len(sys.argv) > 5:
    basis = (float)(sys.argv[5])
  
  print("Parameters:\n"
    "  output_filename: \"{}\"\n"
    "  n_motor_units:    {}\n"
    "  mode:             0 (random placement of fibers)\n"
    "  n_fibers_x:       {} (number of fibers in one coordinate direction) \n"
    "  basis:            {} i.e. exponential distribution will be MU size = {}^mu_no\n".
  format(output_filename, n_motor_units, n_fibers_x, basis, basis))

  n_fibers_y = n_fibers_x
  n_fibers = n_fibers_x * n_fibers_y
  mu_positions = None
  
  # output start and end values of cdf, which should be 0 and 1
  x_list = np.linspace(1, n_motor_units, 100)
  print("cdf(1): {}, cdf(max): {}".format(cdf(1), cdf(n_motor_units)))
  
  # write file
  with open(output_filename,"w") as f:
    values = [draw_sample() for i in range(n_fibers)]
    str_values = map(str, values)
    f.write(" ".join(str_values))
  
# mode 1: centralized placement of fibers around center of MU, assigning all fibers
elif mode == 1:
  
  n_fibers_x = -1
  if len(sys.argv) > 4:
    n_fibers_x = (int)(sys.argv[4])
  else:
    print("The number of fibers in x and y direction has to be known:")
    n_fibers_x = (int)(input("Please enter n_fibers_x and press Enter: "))
  
  if len(sys.argv) > 5:
    basis = (float)(sys.argv[5])
  
  n_max_iterations = 100
  if len(sys.argv) > 6:
    n_max_iterations = (int)(sys.argv[6])
  
  chunk_size = 5          # sizes of the optimization problems that will be solved. A smaller chunk_size means faster optimization problems but more of them
  if len(sys.argv) > 7:
    chunk_size = (int)(sys.argv[7])
  
  print("Parameters:\n"
    "  output_filename: \"{}\"\n"
    "  n_motor_units:    {}\n"
    "  mode:             1 (centralized placement of fibers around center of MU, assigning all fibers)\n"
    "  n_fibers_x:       {} (number of fibers in one coordinate direction) \n"
    "  basis:            {} i.e. exponential distribution will be MU size = {}^mu_no\n"
    "  n_max_iterations: {} (maximum number of iterations of the nonlinear solver, decrease this value if the script takes too long, -1=use pysgpp)\n"
    "  chunk_size:       {} (size of the optimization problems to solve, decrease to get smaller optimization problems but more of them)\n"
    "                       There will be {} optimization problems to be solved\n".
  format(output_filename, n_motor_units, n_fibers_x, basis, basis, n_max_iterations, chunk_size, (int)(np.ceil(n_motor_units / chunk_size))))

  tolerance = 1e-7
  enable_plots = True
  use_accurate_radial_basis_function = False      # True: use Gaussian, this leads to better results, False: use inverse quadratic rbf, this is faster

  # shape parameters for radial basis function
  sigma = 0.1*n_fibers_x    # std deviation for size of motor unit territories
  a_factor = 2/(sigma**2)   # use a higher factor for a smaller width of the RBF, which leads to sharper motor unit territories
  
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

  print("  Randomly chosen MU positions:")  

  # plot actual center points of MUs
  if enable_plots:
    # plot MU centers in 2D plot
    colors = cm.rainbow(np.linspace(0, 1, n_motor_units))

    plt.figure(0)
    for mu_no in range(n_motor_units):
      mu_position = mu_positions[mu_no]
      print("  MU {}, center position: {}".format(mu_no,mu_position))
      x = mu_position[0]
      y = mu_position[1]
    
      color = colors[mu_no,:]
      plt.plot(x,y, '+', markersize=24,color=color)
    filename = "plots/"+output_filename+"_mu_positions.pdf"
    plt.savefig(filename)
    print("  Saved MU positions to file \"{}\".".format(filename))
  
  def pdf_distance(equalization_factors,scaling_factor,i,j,mu_no):
    """ the probability distribution of the motor unit of fiber (i,j) """
    
    # scaling factor would be computed as follows, however, the value is computed beforehand and cached, this is a lot faster
    #scaling_factor = sum([pdf_distance_unscaled(equalization_factors,i,j,t) for t in range(1, n_motor_units+1)])
    
    return pdf_distance_unscaled(equalization_factors,i,j,mu_no) / scaling_factor
    
  def cdf_distance(equalization_factors,i,j,x):
    """ the cdf of pdf_distance """
    
    # compute temporary variable
    scaling_factor = sum([pdf_distance_unscaled(equalization_factors,i,j,t) for t in range(1, n_motor_units+1)])
    
    # discrete case
    x = (int)(np.round(x))
    value = 0
    for t in range(1,x+1):
      value += pdf_distance(equalization_factors,scaling_factor,i,j,t)
    return value
    
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

  def objective(equalization_factors_to_optimize):
    """ objective function for the optimization of the equalization_factors, how well the total fiber distribution follows the exponential pdf() distribution """
    
    # prepare the equalization_factors, part of them are prescribed, part of them are optimized for (equalization_factors_to_optimize)
    for i,index in enumerate(transfer_indices):
      equalization_factors_up_to_chunk[index] = equalization_factors_to_optimize[i]
    
    equalization_factors = equalization_factors_up_to_chunk
    
    # compute scaling factors for given equalization_factors, these scaling factors are only computed once here, instead of every time in pdf_distance
    scaling_factors = np.zeros((n_fibers_y,n_fibers_x))
    for j in range(n_fibers_y):
      for i in range(n_fibers_x):
        scaling_factors[j,i] = sum([pdf_distance_unscaled(equalization_factors,i,j,t) for t in range(1, n_motor_units+1)])
    
    result = 0
    for mu_no in range(1,n_motor_units+1):
      mu_no_total = motor_unit_indices_up_to_chunk[mu_no-1]+1
      #print("mu_no: {}, motor_unit_indices_up_to_chunk: {} -> {}".format(mu_no, motor_unit_indices_up_to_chunk, motor_unit_indices_up_to_chunk[mu_no-1]))
      
      #mu_no_total = mu_no
      required_propability_per_fiber = pdf(mu_no_total)
      
      expected_value_mu = sum([pdf_distance(equalization_factors,scaling_factors[j,i],i,j,mu_no) for j in range(n_fibers_y) for i in range(n_fibers_x)])
      propability_per_fiber = expected_value_mu / n_fibers_total
    
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

  # compute scaling factors for given equalization_factors
  scaling_factors = np.zeros((n_fibers_y,n_fibers_x))
  for j in range(n_fibers_y):
    for i in range(n_fibers_x):
      scaling_factors[j,i] = sum([pdf_distance_unscaled(equalization_factors,i,j,t) for t in range(1, n_motor_units+1)])
  
  # print expected values of number of fiber per MU, before applying the equalization_factors
  print("\nInitial number of fibers for the MUs, this value is not yet exponentially distributed.\n(If you set n_max_iterations very low (which is fast), you'll get MU sizes like this.)")
  n_fibers_total = n_fibers_x * n_fibers_y
  for mu_no in range(1,n_motor_units+1):
    expected_value_mu = sum([pdf_distance(equalization_factors,scaling_factors[j,i],i,j,mu_no) for j in range(n_fibers_y) for i in range(n_fibers_x)])
    print("  MU {}, determined initial number of fibers (expected value): {:.3f}".format(mu_no, expected_value_mu))
    
  # plot radial basis function
  if enable_plots:
    plt.figure(1)
    
    # plot radial basis function
    x_list = np.linspace(mu_positions[0]-n_fibers_x//2, mu_positions[0]+n_fibers_x//2)
    
    previous = use_accurate_radial_basis_function
    
    # plot accurate rbf
    use_accurate_radial_basis_function = True
    y_list = [pdf_distance(unity_equalization_factors,sum([pdf_distance_unscaled(unity_equalization_factors,i,mu_positions[0][1],t) for t in range(1, n_motor_units+1)]),i,mu_positions[0][1],1) for i in x_list]
    plt.plot(x_list, y_list, "b-", label="PDF (Gaussian)")
    y_list = [pdf_distance_unscaled(unity_equalization_factors,i,mu_positions[0][1],1) for i in x_list]
    plt.plot(x_list, y_list, "g-", label="RBF (Gaussian)")
    
    # plot approximated rbf
    use_accurate_radial_basis_function = False
    y_list = [pdf_distance(unity_equalization_factors,sum([pdf_distance_unscaled(unity_equalization_factors,i,mu_positions[0][1],t) for t in range(1, n_motor_units+1)]),i,mu_positions[0][1],1) for i in x_list]
    plt.plot(x_list, y_list, "b--", label="PDF (inverse quadratic)")
    y_list = [pdf_distance_unscaled(unity_equalization_factors,i,mu_positions[0][1],1) for i in x_list]
    plt.plot(x_list, y_list, "g--", label="RBF (inverse quadratic)")
    use_accurate_radial_basis_function = previous
    
    plt.title("Comparison of Gaussian and inverse quadratic RBFs")
    plt.legend(loc="best")
    
  # find equalization_factors that fulfill the required exponential distribution of MU sizes
  print("\nNow run optimizer. Optimize factors to obtain exponential distribution of MU sizes. Maximum number of iterations: {}".format(n_max_iterations))
    
  def optimize():
    """
    solve the optimization problem for the equalization_factors, for given motor unit positions
    """
    
    equalization_factors_to_optimize = [1 for _ in range(n_motor_units_chunk)]
    print("optimize(n_motor_units={}, len(mu_positions)={}, n_motor_units_chunk={})".format(n_motor_units, len(mu_positions), n_motor_units_chunk))
  
    # solve optimization problem using sparse grids, if this fails, fall back to scipy.optimize
    try:
      import pysgpp
      import traceback
      print("The module pysgpp is available.")
      
      if n_max_iterations != -1:
        raise Exception("The value of <n_max_iterations> is not set to -1, this means do not use pysgpp but scipy.optimize.")

      class ExampleFunction(pysgpp.ScalarFunction):
        """objective function"""
        def __init__(self):
          super(ExampleFunction, self).__init__(size_optimization)
        
        def eval(self, x):
          """Evaluates the function."""
          print("x:{}".format(x))
          try:
            #scaled_x = [value*2 for value in x]
            for i in range(len(x)):
              x[i] = x[i]*2
            result = objective(x)
          except Exception as a:
            print(a)
            traceback.print_exc() 
          print("result: {}".format(result))
          return result
        
      #pysgpp.omp_set_num_threads(64)
      pysgpp.omp_set_num_threads(2)
      
      # increase verbosity of the output
      pysgpp.Printer.getInstance().setVerbosity(2)
      
      # objective function
      f = ExampleFunction()
      
      # dimension of domain
      d = f.getNumberOfParameters()
      
      # B-spline degree
      p = 3
      
      # maximal number of grid points
      N = 30
      
      # adaptivity of grid generation
      gamma = 0.95
      
      ## First, we define a grid with modified B-spline basis functions and
      ## an iterative grid generator, which can generate the grid adaptively.
      grid = pysgpp.Grid.createModBsplineGrid(d, p)
      gridGen = pysgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, gamma)

      ## With the iterative grid generator, we generate adaptively a sparse grid.
      print("Generating grid...\n")

      if not gridGen.generate():
        print("Grid generation failed, exiting.")
        sys.exit(1)

      ## Then, we hierarchize the function values to get hierarchical B-spline
      ## coefficients of the B-spline sparse grid interpolant
      ## \f$\tilde{f}\colon [0, 1]^d \to \mathbb{R}\f$.
      print("Hierarchizing...\n")
      functionValues = gridGen.getFunctionValues()
      coeffs = pysgpp.DataVector(len(functionValues))
      hierSLE = pysgpp.HierarchisationSLE(grid)
      sleSolver = pysgpp.AutoSLESolver()

      # solve linear system
      if not sleSolver.solve(hierSLE, gridGen.getFunctionValues(), coeffs):
        print("Solving failed, exiting.")
        sys.exit(1)

      ## We define the interpolant \f$\tilde{f}\f$ and its gradient
      ## \f$\nabla\tilde{f}\f$ for use with the gradient method (steepest descent).
      ## Of course, one can also use other optimization algorithms from
      ## sgpp::optimization::optimizer.
      print("Optimizing smooth interpolant...\n")
      ft = pysgpp.InterpolantScalarFunction(grid, coeffs)
      ftGradient = pysgpp.InterpolantScalarFunctionGradient(grid, coeffs)
      gradientDescent = pysgpp.OptGradientDescent(ft, ftGradient)
      x0 = pysgpp.DataVector(d)

      ## The gradient method needs a starting point.
      ## We use a point of our adaptively generated sparse grid as starting point.
      ## More specifically, we use the point with the smallest
      ## (most promising) function value and save it in x0.
      gridStorage = gridGen.getGrid().getStorage()

      # index of grid point with minimal function value
      x0Index = 0
      fX0 = functionValues[0]
      for i in range(1, len(functionValues)):
        if functionValues[i] < fX0:
          fX0 = functionValues[i]
          x0Index = i

      x0 = gridStorage.getCoordinates(gridStorage.getPoint(x0Index));
      ftX0 = ft.eval(x0)

      print("x0 = {}".format(x0))
      print("f(x0) = {:.6g}, ft(x0) = {:.6g}\n".format(fX0, ftX0))

      ## We apply the gradient method and print the results.
      gradientDescent.setStartingPoint(x0)
      gradientDescent.optimize()
      xOpt = gradientDescent.getOptimalPoint()
      ftXOpt = gradientDescent.getOptimalValue()
      fXOpt = f.eval(xOpt)

      print("\nxOpt = {}".format(xOpt))
      print("f(xOpt) = {:.6g}, ft(xOpt) = {:.6g}\n".format(fXOpt, ftXOpt))

      scaled_x = [value*2 for value in xOpt]
      equalization_factors_to_optimize = scaled_x
      
    except Exception as e:
      print(e)
      
      bounds = [(0,None) for i in range(n_motor_units_chunk)]

      result = scipy.optimize.minimize(objective, equalization_factors_to_optimize, bounds=bounds, options={"disp": True, "maxiter": n_max_iterations})
      print("Result of optimizer: success: {}, status: {}, message: {}, n objective evaluations: {}, n iterations: {}".format(result.success, result.status, result.message, result.nfev, result.nit))
      equalization_factors_to_optimize = result.x
        
    return equalization_factors_to_optimize
  
  # divide the problem in chunks of 10 motor units that will be optimized for
  n_chunks = (int)(np.ceil(n_motor_units / chunk_size))
  
  n_motor_units_total = n_motor_units
  mu_positions_total = mu_positions.copy()
  
  equalization_factors_total = [1 for i in range(n_motor_units)]
  motor_unit_indices_collected = []
  
  # set initial values
  equalization_factors = [1 for i in range(n_motor_units)]
  n_motor_units_up_to_chunk = 0
  motor_unit_indices_up_to_chunk = []
  
  for chunk_no in range(n_chunks):
    # select motor units
    
    motor_unit_indices_chunk = list(range(chunk_no,n_motor_units_total,n_chunks))     # indices into list of all motor units
    motor_unit_indices_up_to_chunk += motor_unit_indices_chunk                        # indices into list of all motor units
    motor_unit_indices_up_to_chunk = sorted(motor_unit_indices_up_to_chunk)
    
    transfer_indices = []
    for index in motor_unit_indices_chunk:
      i = motor_unit_indices_up_to_chunk.index(index)
      transfer_indices.append(i)
    
    n_motor_units_chunk = len(motor_unit_indices_chunk)
    n_motor_units_up_to_chunk += n_motor_units_chunk
    
    mu_positions_up_to_chunk = [mu_positions_total[index] for index in motor_unit_indices_up_to_chunk]
    
    # compute scaling factor for the pdf
    scaling_factor_pdf = 0
    for mu_no in range(1,n_motor_units_up_to_chunk+1):
      mu_no_total = motor_unit_indices_up_to_chunk[mu_no-1]+1
      scaling_factor_pdf += pdf_unscaled(mu_no_total)
    
    # prepare equalization factors
    equalization_factors_up_to_chunk = [equalization_factors_total[index] for index in motor_unit_indices_up_to_chunk]
    
    print("chunk #{}/{}, n MUs chunk: {}, so far: {}, motor_unit_indices_chunk: {}, motor_unit_indices_up_to_chunk: {}, transfer_indices: {}".\
    format(chunk_no, n_chunks, n_motor_units_chunk, n_motor_units_up_to_chunk, motor_unit_indices_chunk, motor_unit_indices_up_to_chunk, transfer_indices))
  
    # run optimization
    # needs:
    # equalization_factors_to_optimize      (initialized optimization vector)
    # equalization_factors_up_to_chunk      (all equalization foctors up to the chunk)
    # n_motor_units        (number of motor units up to chunk)
    # mu_positions         (corresponding mu positions)
    # motor_unit_transfer_indices   (indices into equalization_factors_up_to_chunk of the currently considered motor units)
    # scaling_factor_pdf = sum([pdf_unscaled(x) for x in range(1,n_motor_units+1)])
    
    n_motor_units = n_motor_units_up_to_chunk
    mu_positions = mu_positions_up_to_chunk
   
    equalization_factors_chunk = optimize()
    print("equalization_factors_chunk: {}".format(equalization_factors_chunk))
    
    # store the newly found equalization_factors_chunk to equalization_factors_total
    for i,index in enumerate(motor_unit_indices_chunk):
      equalization_factors_total[index] = equalization_factors_chunk[i]
    
    print("  equalization_factors_total: {}".format(equalization_factors_total))
  
    # the rest in this loop is only for debugging output
    if True:
      # prepare the equalization_factors, part of them are prescribed, part of them are optimized for (equalization_factors_chunk)
      for i,index in enumerate(transfer_indices):
        equalization_factors_up_to_chunk[index] = equalization_factors_chunk[i]
      
      # compute scaling factors for given equalization_factors, these scaling factors are only computed once here, instead of every time in pdf_distance
      scaling_factors = np.zeros((n_fibers_y,n_fibers_x))
      for j in range(n_fibers_y):
        for i in range(n_fibers_x):
          scaling_factors[j,i] = sum([pdf_distance_unscaled(equalization_factors_up_to_chunk,i,j,t) for t in range(1, n_motor_units+1)])
      
      result = 0
    
      # print the expected and actual number of fibers per MU
      print("\n  MU sizes, expected values from probability and actual realized sizes:")
      sum_optimal_p = 0
      sum_determined_p = 0
      for mu_no in range(1,n_motor_units+1):
        expected_value_mu = sum([pdf_distance(equalization_factors_up_to_chunk,scaling_factors[j,i],i,j,mu_no) for j in range(n_fibers_y) for i in range(n_fibers_x)])
        
        mu_no_total = motor_unit_indices_up_to_chunk[mu_no-1]+1
        required_propability_per_fiber = pdf(mu_no_total)
        #required_propability_per_fiber = 1./n_motor_units
        
        expected_value_mu = sum([pdf_distance(equalization_factors_up_to_chunk,scaling_factors[j,i],i,j,mu_no) for j in range(n_fibers_y) for i in range(n_fibers_x)])
        propability_per_fiber = expected_value_mu / n_fibers_total
      
        penalty = (required_propability_per_fiber - propability_per_fiber)**2
        result += penalty
        sum_optimal_p += pdf(mu_no_total)*n_fibers_total
        sum_determined_p += expected_value_mu
        print("  MU {}, n_fibers by optimal p: {:.3f}, by determined p: {:.3f}, penalty: {:.3e}=({:.3e}-{:.3e})^2".
          format(mu_no, pdf(mu_no_total)*n_fibers_total, expected_value_mu, penalty, required_propability_per_fiber, propability_per_fiber))
    
      print("  final objective: {}".format(result))
      print("  sum:    {}, {}, (total number of fibers: {})\n".format(sum_optimal_p, sum_determined_p, n_fibers_total))
  
  n_motor_units = n_motor_units_total
  mu_positions = mu_positions_total
  equalization_factors = equalization_factors_total
  scaling_factor_pdf = sum([pdf_unscaled(x) for x in range(1,n_motor_units+1)])
  
  # compute scaling factors for given equalization_factors
  scaling_factors = np.zeros((n_fibers_y,n_fibers_x))
  for j in range(n_fibers_y):
    for i in range(n_fibers_x):
      scaling_factors[j,i] = sum([pdf_distance_unscaled(equalization_factors,i,j,t) for t in range(1, n_motor_units+1)])
  
  t_end = time.time()
  print("\nTotal duration: {:.2f} s\n\n".format(t_end-t_start))
  
  # print the result
  #print("determined equalization_factors:")
  #print(equalization_factors)
    
  # initialize i,j to the center fiber in the muscle
  i = (int)(n_fibers_x/2)
  j = (int)(n_fibers_y/2)
  
  n_points_to_plot = 100
  
  # output start and end values of cdf, which should be 0 and 1
  x_list = list(range(1,n_motor_units+1)) #np.linspace(1,n_motor_units,n_points_to_plot)
  #print("at center, (i,j) = ({},{}): cdf_distance(0): {} (should be 0), cdf_distance(max): {} (should be 1)".
  #  format(i,j,cdf_distance(equalization_factors,i,j,0), cdf_distance(equalization_factors,i,j,n_motor_units)))

  if enable_plots:
    
    # plot pdf
    plt.figure(2)
    # compute scaling factors for unity equalization_factors
    unity_scaling_factors = np.zeros((n_fibers_y,n_fibers_x))
    for j in range(n_fibers_y):
      for i in range(n_fibers_x):
        unity_scaling_factors[j,i] = sum([pdf_distance_unscaled(unity_equalization_factors,i,j,t) for t in range(1, n_motor_units+1)])
    
    #print("test if inverse is correct: 0=",max([inverse_cdf_distance(i,j,cdf_distance(i,j,x))-x for x in x_list]))
    plt.plot(x_list, [pdf_distance(unity_equalization_factors,unity_scaling_factors[j,i],i,j,x) for x in x_list], '-o', label="pdf")
    plt.gca().set_ylim(0,1)
    plt.xticks(range(n_motor_units+1))
    plt.ylabel("probability for a fiber to be in this MU\n (only according to distance)")
    plt.grid(which='both')
    plt.legend()
    plt.savefig("plots/"+output_filename+"_pdf_distance.pdf")

    # plot pdf total
    plt.figure(3)
    #print("test if inverse is correct: 0=",max([inverse_cdf_distance(i,j,cdf_distance(i,j,x))-x for x in x_list]))
    plt.plot(x_list, [pdf_distance(equalization_factors,scaling_factors[j,i],i,j,x) for x in x_list], '-o', label="pdf")
    plt.gca().set_ylim(0,1)
    plt.xticks(range(n_motor_units+1))
    plt.ylabel("probability for a fiber to be in this MU")
    plt.grid(which='both')
    plt.legend()
    plt.savefig("plots/"+output_filename+"_pdf_total.pdf")

    # plot cdf and inverse
    plt.figure(4)
    plt.plot(x_list, [cdf_distance(equalization_factors,i,j,x) for x in x_list], '-', label="cdf")
    y_list = np.linspace(0,1,n_points_to_plot)
    plt.plot(y_list, [inverse_cdf_distance(equalization_factors,i,j,y) for y in y_list], '-o', label="inverse cdf")
    plt.xlabel("motor unit no")
    plt.xticks(range(n_motor_units+1))
    plt.gca().axis('equal')
    plt.grid(which='both')
    #plt.plot(x_list, [inverse_distribution(x) for x in x_list], label="f^-1")
    plt.legend()
    plt.savefig("plots/"+output_filename+"_cdf.pdf")
    #plt.show()
  
  # create actual motor unit numbers by sampling according to the pdf
  mu_nos_for_fibers = []
  
  # try sampling the probability 10 times and take the best match
  print("Sampling MUs for fibers, best out of 5")
  best_error = None
  for iteration_no in range(5):
    
    # loop over fibers
    iteration_mu_nos_for_fibers = []
    for j in range(n_fibers_y):
      for i in range(n_fibers_x):
        mu_no = draw_sample(i,j)
        if mu_no > n_motor_units:
          mu_no = n_motor_units
        iteration_mu_nos_for_fibers.append(mu_no)
    
    bins = [x-0.5 for x in range(1,n_motor_units+2)]
    numbers,edges = np.histogram(iteration_mu_nos_for_fibers,bins)
    
    # assert that there are no motor units without any fibers
    for mu_no in range(1,n_motor_units+1):
      if numbers[mu_no-1] == 0:
        
        # determine probabilities
        mu_probabilities = []
        for j in range(n_fibers_y):
          for i in range(n_fibers_x):
            # determine probability for MU mu_no at fiber (i,j)
            scaling_factor = sum([pdf_distance_unscaled(equalization_factors,i,j,t) for t in range(1, n_motor_units+1)])
            p = pdf_distance(equalization_factors,scaling_factor,i,j,mu_no)
            mu_probabilities.append(p)
            
        max_p = max(mu_probabilities)
        index = mu_probabilities.index(max_p)
        
        # set this fiber to be of the MU
        iteration_mu_nos_for_fibers[index] = mu_no
        
        # create new histogram
        numbers,edges = np.histogram(iteration_mu_nos_for_fibers,bins)
  
    error = 0
    for mu_no in range(1,n_motor_units+1):
      error += (numbers[mu_no-1] - pdf(mu_no)*n_fibers_total)**2
    print("  error: {:.1f} (number fibers per MU: {})".format(error, numbers))
    if best_error is None or error < best_error:
      best_error = error
      mu_nos_for_fibers = list(iteration_mu_nos_for_fibers)
      
  bins = [x-0.5 for x in range(1,n_motor_units+2)]
  numbers,edges = np.histogram(mu_nos_for_fibers,bins)

  # print the expected and actual number of fibers per MU
  print("\nMU sizes, number of fibers expected by optimal probability, number of fibers expected by determined probability and actual realized sizes:")
  sum_optimal_p = 0
  sum_determined_p = 0
  for mu_no in range(1,n_motor_units+1):
    expected_value_mu = sum([pdf_distance(equalization_factors,scaling_factors[j,i],i,j,mu_no) for j in range(n_fibers_y) for i in range(n_fibers_x)])
    sum_optimal_p += pdf(mu_no)*n_fibers_total
    sum_determined_p += expected_value_mu
    print("  MU {}, n. fibers optimal: {:.3f}, n. fibers expected: {:.3f}, realized: {}".
      format(mu_no, pdf(mu_no)*n_fibers_total, expected_value_mu, numbers[mu_no-1]))
  print("sum:    {}, {}, (total number of fibers: {})".format(sum_optimal_p, sum_determined_p, n_fibers_total))

      
  # write all motor units to output file
  with open(output_filename,"w") as f:
    
    str_values = map(str, mu_nos_for_fibers)
    f.write(" ".join(str_values))
  
# mode 2: centralized placement of fibers around center of MU, not assigning all fibers, this is analogous to multidomain
elif mode == 2:

  n_fibers_x = -1
  if len(sys.argv) > 4:
    n_fibers_x = (int)(sys.argv[4])
  else:
    print("The number of fibers in x and y direction has to be known:")
    n_fibers_x = (int)(input("Please enter n_fibers_x and press Enter: "))
  
  if len(sys.argv) > 5:
    basis = (float)(sys.argv[5])
  
  print("Parameters:\n"
    "  output_filename: \"{}\"\n"
    "  n_motor_units:    {}\n"
    "  mode:             2 (centralized placement of fibers around center of MU, not assigning all fibers)\n"
    "  n_fibers_x:       {} (number of fibers in one coordinate direction) \n"
    "  basis:            {} i.e. exponential distribution will be MU size = {}^mu_no\n".
  format(output_filename, n_motor_units, n_fibers_x, basis, basis))

  enable_plots = True
  use_accurate_radial_basis_function = False      # True: use Gaussian, this leads to better results, False: use inverse quadratic rbf, this is faster

  # shape parameters for radial basis function
  sigma = 0.1*n_fibers_x
  a_factor = 2/(sigma**2)   # use a higher factor for a smaller width of the RBF, which leads to sharper motor unit territories
  
  n_fibers_y = n_fibers_x
  n_fibers_total = n_fibers_x * n_fibers_y

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

  print("  Randomly chosen MU positions:")  

  # plot actual center points of MUs
  if enable_plots:
    # plot MU centers in 2D plot
    colors = cm.rainbow(np.linspace(0, 1, n_motor_units))

    plt.figure(0)
    for mu_no in range(n_motor_units):
      mu_position = mu_positions[mu_no]
      print("  MU {}, center position: {}".format(mu_no,mu_position))
      x = mu_position[0]
      y = mu_position[1]
    
      color = colors[mu_no,:]
      plt.plot(x,y, '+', markersize=24,color=color)
    filename = "plots/"+output_filename+"_mu_positions.pdf"
    plt.savefig(filename)
    print("  Saved MU positions to file \"{}\".".format(filename))
  
  equalization_factors = [1 for _ in range(n_motor_units)]
  
  mu_nos_for_fibers = []
  
  def pdf_total(i,j,mu_no):
    """ probability for fiber (i,j) to be of motor unit mu_no """
    return pdf_distance_unscaled(equalization_factors,i,j,mu_no) * pdf(mu_no)
  
  # at first check maximum probability value that would occur in the sum
  p_scaling = max([sum([pdf_total(i,j,mu_no) for mu_no in range(1,n_motor_units+1)]) for j in range(n_fibers_y) for i in range(n_fibers_x)])
  
  print("p_scaling: {}".format(p_scaling))
  print("The fibers that are not assigned to any MU get MU no. {}.".format(n_motor_units+1))
  
  # scale down all probabilities by this value
  
  print("Sampling MUs for fibers, best out of 20")
  best_error = None
  for iteration_no in range(1):
      
    iteration_mu_nos_for_fibers = []
    
    # loop over fibers
    for j in range(n_fibers_y):
      for i in range(n_fibers_x):
        
        # loop over motor units and get probabilities
        mu_probabilities = []
        for mu_no in range(1,n_motor_units+1):
          p = pdf_total(i,j,mu_no) / p_scaling
          mu_probabilities.append(p)
          
        missing_p = 1 - sum(mu_probabilities)
        mu_probabilities.append(missing_p)
          
        # draw motor unit
        mu_no = np.random.choice(list(range(1,n_motor_units+2)), p=mu_probabilities)
        iteration_mu_nos_for_fibers.append(mu_no)
   
    bins = [x-0.5 for x in range(1,n_motor_units+2)]
    numbers,edges = np.histogram(iteration_mu_nos_for_fibers,bins)
    
    # assert that there are no motor units without any fibers
    for mu_no in range(1,n_motor_units+1):
      if numbers[mu_no-1] == 0:
        
        # determine probabilities
        mu_probabilities = []
        for j in range(n_fibers_y):
          for i in range(n_fibers_x):
            # determine probability for MU mu_no at fiber (i,j)
            p = pdf_total(i,j,mu_no) / p_scaling
            mu_probabilities.append(p)
            
        max_p = max(mu_probabilities)
        index = mu_probabilities.index(max_p)
        
        # set this fiber to be of the MU
        iteration_mu_nos_for_fibers[index] = mu_no
        
        # create new histogram
        numbers,edges = np.histogram(iteration_mu_nos_for_fibers,bins)
  
    error = 0
    for mu_no in range(1,n_motor_units+1):
      error += (numbers[mu_no-1] - pdf(mu_no)*n_fibers_total)**2
    #print("  error: {:.1f} (number fibers per MU: {})".format(error, numbers))
    if best_error is None or error < best_error:
      best_error = error
      mu_nos_for_fibers = list(iteration_mu_nos_for_fibers)
      
  # determine number of unassigned fibers
  n_unassigned_fibers = sum([1 if mu_nos_for_fibers[i] == n_motor_units+1 else 0 for i in range(len(mu_nos_for_fibers))])
      
  print("{} of {} fibers, i.e. {:.0f} % have no motor unit.".format(n_unassigned_fibers, n_fibers_total, 100*n_unassigned_fibers/n_fibers_total))
  
  bins = [x-0.5 for x in range(1,n_motor_units+2)]
  numbers,edges = np.histogram(mu_nos_for_fibers,bins)

  # print the expected and actual number of fibers per MU
  print("\nMU sizes:")
  for mu_no in range(1,n_motor_units+1):
    print("  MU {}: {} fiber(s)".
      format(mu_no, numbers[mu_no-1]))

  # write all motor units to output file
  with open(output_filename,"w") as f:
    
    str_values = map(str, mu_nos_for_fibers)
    f.write(" ".join(str_values))

# plot results, independent of mode
# ---------------------------------

# parse generated mu assignments
f = open(output_filename,"r")
line = f.readline()
mu_nos_for_fibers = line.split(" ")
mu_nos_for_fibers = list(map(int, mu_nos_for_fibers))
  
# plot fibers in 2D plot
#print("n_motor_units: {}, min: {}".format(n_motor_units, (int)(min(mu_nos_for_fibers))))

colors = cm.rainbow(np.linspace(0, 1, n_motor_units))
color_invalid_mu = np.reshape(np.array((1,1,1,1)), (1,4))
colors = np.concatenate((colors, color_invalid_mu))

index = 0
point_colors = []
for j in range(n_fibers_y):
  for i in range(n_fibers_x):
    mu_no = (int)(mu_nos_for_fibers[index])
    point_colors.append(colors[mu_no-1,:])
    #plt.plot(i,j,'o',color=colors[mu_no-1,:])
    
    index += 1
    
fig = plt.figure(5,figsize=(8,8))
  
# plot actual center points of MUs
if mu_positions is not None:
  for mu_no in range(n_motor_units):
    mu_position = mu_positions[mu_no]
    x = mu_position[0]
    y = mu_position[1]
    
    color = colors[mu_no,:]
    #print(mu_position,x,y,color)
    plt.plot(x,y, 'x', markersize=24,markeredgewidth=2,color=color)

X,Y = np.meshgrid(range(n_fibers_x),range(n_fibers_y))
plt.scatter(X,Y,color=point_colors,marker="s")
           
ax = plt.gca()

legend_elements = []

for mu_no in range(n_motor_units):
  legend_elements.append(Line2D([0], [0], marker="s", color=colors[mu_no], lw=0, label='MU {}'.format(mu_no+1)))
        
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
plt.axis('equal')
plt.tight_layout()
plt.savefig("plots/"+output_filename+"_2d_fiber_distribution.pdf")

# plot total distribution of MUs

fig = plt.figure(6)
bins = [x-0.5 for x in range(1,n_motor_units+2)]
numbers,edges = np.histogram(mu_nos_for_fibers,bins)

plt.hist(mu_nos_for_fibers,bins=bins, align='mid', rwidth=0.8)
a = numbers[-1] / (basis**n_motor_units)
xlist = np.linspace(1,n_motor_units,5*n_motor_units)
plt.plot(xlist, [basis**x * a for x in xlist], label='${}^x$'.format(basis))
#plt.plot(xlist, [np.exp(x)/np.exp(n_motor_units)*bins[-1] for x in xlist])
ax = fig.gca()
ax.set_xlim(1,n_motor_units)
from matplotlib.ticker import MaxNLocator
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.xlabel("motor unit index")
plt.ylabel("count")
plt.legend()
plt.savefig("plots/"+output_filename+"_fiber_distribution.pdf")

print("\nResult plots were written in plot subdirectory as \"plots/"+output_filename+"_*.pdf\".")
#plt.show()

