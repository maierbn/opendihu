# Diffusion 1D1
#
# command arguments: <name> <number elements>

import numpy as np
#import scipy.integrate
import sys

n = 5   # number of elements
nx = 4*n
ny = 3*n
name = ""

if len(sys.argv) > 1:
  if "check_results.py" not in sys.argv[0]:
    name = sys.argv[0]
    n = int(sys.argv[1])
    
    print("name: \"{}\", n: {}".format(name, n))

def initial_values_function(x,y):
  if 0.75 < x < 1.25 and 0.75 < y < 1.75:
    return np.sin((x-0.75)/0.5*np.pi) * np.sin((y-0.75)/1.0*np.pi)
  if 1.5 < x < 3.0 and 1.0 < y < 2.25:
    return np.sin((x-1.5)/1.5*np.pi) * np.sin((y-1.0)/1.25*np.pi)*0.5
  return 0
  
end_time = 0.5
c = 0.2     # prefactor in equation

if "hermite" in name:
  end_time = 0.1

# set node positions for unstructured mesh
n_nodes_x = nx+1
n_nodes_y = ny+1
if "quadratic" in name:
  n_nodes_x = 2*nx + 1
  n_nodes_y = 2*ny + 1
node_positions = [[x,y] for y in np.linspace(0,3.0,n_nodes_y) for x in np.linspace(0,4.0,n_nodes_x) ]
hx = 4.0/(n_nodes_x)
hy = 4.0/(n_nodes_y)

# move node positions for unstructured grid
if "unstructured" in name:
  for iy in range(n_nodes_y):
    for ix in range(n_nodes_x):
      i = iy*n_nodes_x + ix
      node_positions[i][0] += hx*np.sin(float(ix)/n_nodes_x*4*np.pi)*0.1
      node_positions[i][1] += hy*np.sin(float(iy)/n_nodes_y*6*np.pi)*0.1

# set elements for unstructured mesh
elements = []
for iy in range(n_nodes_y-1):
  for ix in range(n_nodes_x-1):
    i = iy*n_nodes_x + ix
    elements.append([i,i+1,i+n_nodes_x,i+n_nodes_x+1])

# elements for quadratic unstructured
if name == "2d_unstructured_deformable_quadratic":
  elements = []
  for iy in range(n_nodes_y-1,2):
    for ix in range(n_nodes_x-1,2):
      i = iy*n_nodes_x + ix
      elements.append([i,i+1,i+2, i+n_nodes_x,i+n_nodes_x+1,i+n_nodes_x+2, i+2*n_nodes_x,i+2*n_nodes_x+1,i+2*n_nodes_x+2])

# set initial values
initial_values = [initial_values_function(x0,x1) for x1 in np.linspace(0,3.0,n_nodes_y) for x0 in np.linspace(0,4.0,n_nodes_x) ]

#print "initial_values: ",initial_values
#print "len: ", len(initial_values)


config = {
  "ExplicitEuler" : {
    "initialValues": initial_values,
    "numberTimeSteps": 1500,
    "endTime": end_time,
    "FiniteElementMethod" : {
      "nElements": [nx,ny],
      "nodeDimension": 1,
      "nodePositions": node_positions,
      "elements": elements,
      "physicalExtent": [4.0,3.0],
      "relativeTolerance": 1e-15,
      "prefactor": c
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      #{"format": "ExFile", "filename": "out/"+name, "outputInterval": 2},
      {"format": "PythonFile", "filename": "out/"+name, "outputInterval": 5, "binary":False, "onlyNodalValues":True},
    ]
  },
}

# output config in a readable format
if False:
  import pprint 
  pp = pprint.PrettyPrinter()
  pp.pprint(config)
