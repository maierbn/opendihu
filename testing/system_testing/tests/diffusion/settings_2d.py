# Diffusion 2D
#
# command arguments: <name> <number elements>

import numpy as np
#import scipy.integrate
import sys

n = 5   # number of elements
nx = 0
ny = 0
name = ""

if len(sys.argv) > 1:
  if "check_results.py" not in sys.argv[0]:
    name = sys.argv[0]
    n = int(sys.argv[1])
    nx = int(np.sqrt(n))
    ny = nx+1
    
    print("name: \"{}\", n: {}, nx: {}, ny: {}".format(name, n, nx, ny))

# function for initial values
def initial_values_function(x,y):
  if 0.75 < x < 1.25 and 0.75 < y < 1.75:
    return np.sin((x-0.75)/0.5*np.pi) * np.sin((y-0.75)/1.0*np.pi)
  if 1.5 < x < 3.0 and 1.0 < y < 2.25:
    return np.sin((x-1.5)/1.5*np.pi) * np.sin((y-1.0)/1.25*np.pi)*0.5
  return 0
  
end_time = 0.5
c = 0.2     # prefactor in equation

if "hermite" in name:
  end_time = 0.05

n_nodes_x = 0   # this needs to be defined such that check_results.py works
n_nodes_y = 0  
if "linear" in name or "hermite" in name:
  n_nodes_x = nx + 1
  n_nodes_y = ny + 1

if "quadratic" in name:
  n_nodes_x = 2*nx + 1
  n_nodes_y = 2*ny + 1

physicalExtent = [4.0, 3.0]

if "fixed" in name:
  physicalExtent = [nx, ny]

hx = physicalExtent[0]/(n_nodes_x-1.0)
hy = physicalExtent[1]/(n_nodes_y-1.0)

# set node positions for unstructured mesh
node_positions = []
for y in list(np.linspace(0,physicalExtent[1],n_nodes_y)):
  for x in list(np.linspace(0,physicalExtent[0],n_nodes_x)):
    node_positions.append([x,y])
  
# move node positions
move_factor = 0.1
for i in range(len(node_positions)):
  node_positions[i][0] += move_factor * hx * np.sin(float(i%n_nodes_x)/n_nodes_x*4*np.pi)
  node_positions[i][1] += move_factor * hy * np.sin(float(int(i/n_nodes_x))/n_nodes_y*4*np.pi)

# set elements for unstructured mesh
elements = []
if "linear" in name or "hermite" in name:
  for iy in range(n_nodes_y-1):
    for ix in range(n_nodes_x-1):
      i0 = iy*n_nodes_x + ix
      i1 = iy*n_nodes_x + ix+1
      i2 = (iy+1)*n_nodes_x + ix
      i3 = (iy+1)*n_nodes_x + ix+1
      
      elements.append([i0, i1, i2, i3])
      
elif "quadratic" in name:
  for iy in range(n_nodes_y-2):
    for ix in range(n_nodes_x-2):
      i0 = iy*n_nodes_x + ix
      i1 = iy*n_nodes_x + ix+1
      i2 = iy*n_nodes_x + ix+2
      i3 = (iy+1)*n_nodes_x + ix
      i4 = (iy+1)*n_nodes_x + ix+1
      i5 = (iy+1)*n_nodes_x + ix+2
      i6 = (iy+2)*n_nodes_x + ix
      i7 = (iy+2)*n_nodes_x + ix+1
      i8 = (iy+2)*n_nodes_x + ix+2
      
      elements.append([i0, i1, i2, i3, i4, i5, i6, i7, i8])
      
# set initial values
initial_values = [initial_values_function(x0,x1) for x1 in np.linspace(0,physicalExtent[1],n_nodes_y) for x0 in np.linspace(0,physicalExtent[0],n_nodes_x) ]

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
      "physicalExtent": physicalExtent,
      "relativeTolerance": 1e-15,
      "prefactor": c
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      #{"format": "ExFile", "filename": "out/"+name, "outputInterval": 2},
      {"format": "PythonFile", "filename": "out/"+name, "outputInterval": 10, "binary":False, "onlyNodalValues":True},
    ]
  },
}

# output config in a readable format
if False:
  import pprint 
  pp = pprint.PrettyPrinter()
  pp.pprint(config)
