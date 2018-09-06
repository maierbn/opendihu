# Diffusion 1D
#
# command arguments: <name> <number elements>

import numpy as np
#import scipy.integrate
import sys

n = 20   # number of elements
name = ""

if len(sys.argv) > 1:
  if "check_results.py" not in sys.argv[0]:
    name = sys.argv[0]
    n = int(sys.argv[1])
    
    print("name: \"{}\", n: {}".format(name, n))

def initial_values_function(x):
  if 0.25 < x < 0.75:
    return np.sin((2*x-0.5)*np.pi)+0.5
  if 1.5 < x < 3.0:
    return np.sin((x-1.5)/1.5*np.pi)*0.5+0.5
  return 0+0.5
  
end_time = 0.5
c = 0.6     # prefactor in equation

if "hermite" in name:
  end_time = 0.1

# set node positions for unstructured mesh
n_nodes = n+1
if "quadratic" in name:
  n_nodes = 2*n + 1
node_positions = list(np.linspace(0,4.0,n_nodes))
h = 4.0/(n_nodes)

# move node positions
for i in range(n_nodes):
  node_positions[i] += h*np.sin(float(i)/n_nodes*4*np.pi)

# set elements for unstructured mesh
elements = [[i, i+1] for i in range(n_nodes-1)]

if name == "1d_unstructured_deformable_quadratic":
  elements = [[i, i+1, i+2] for i in range(0,n_nodes-1,2)]

# set initial values
initial_values = [initial_values_function(x) for x in np.linspace(0,4.0,n_nodes)]

  

config = {
  "ExplicitEuler" : {
    "initialValues": initial_values,
    "numberTimeSteps": 1500,
    "endTime": end_time,
    "FiniteElementMethod" : {
      "nElements": n,
      "nodeDimension": 1,
      "nodePositions": node_positions,
      "elements": elements,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "prefactor": c
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      {"format": "ExFile", "filename": "out/"+name, "outputInterval": 2},
      {"format": "PythonFile", "filename": "out/"+name, "outputInterval": 2, "binary":False, "onlyNodalValues":True},
    ]
  },
}
