# Diffusion 1D
#
# command arguments: <name> <number elements>

import numpy as np
import scipy.integrate
import sys

n = 20   # number of elements
name = ""

if len(sys.argv) > 1:
  if "check_results.py" not in sys.argv[0]:
    name = sys.argv[0]
    n = int(sys.argv[1])
    
    print "name: \"{}\", n: {}".format(name, n)

def initial_values_function(x):
  if 0.25 < x < 0.75:
    return np.sin((2*x-0.5)*np.pi)
  if 1.5 < x < 3.0:
    return np.sin((x-1.5)/1.5*np.pi)*0.5
  return 0
  
# set initial values
initial_values = [initial_values_function(x) for x in np.linspace(0,4.0,n+1)]

# set node positions for unstructured mesh
node_positions = np.linspace(0,4.0,n+1)
h = 4.0/(n+1)

# move node positions
for i in range(n+1):
  node_positions[i] += h*np.sin(float(i)/(n+1)*4*np.pi)

# define analytic solution for testing
def analytic_solution(x,t):
  
  def h(x,t):
    c = 1.0
    return 1./np.sqrt(4*np.pi*c*t)*np.exp(-x**2 / (4*c*t))
  
  def integrand(y):
    return initial_values_function(y)*h(x-y,t)
    
  (value,error) = scipy.integrate.quad(integrand, 0.0, 4.0)
  return value
  

config = {
  "ExplicitEuler" : {
    "initialValues": initial_values,
    "numberTimeSteps": 1000,
    "endTime": 0.1,
    "FiniteElementMethod" : {
      "nElements": n,
      "nodePositions": node_positions,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/"+name, "outputInterval": 2, "binary":False, "onlyNodalValues":True}
    ]
  },
}
