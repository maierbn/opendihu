# Laplace 2D
#
# command arguments: <name> <number elements>

import numpy as np
#import scipy.integrate
import sys

nx = 5   # number of elements
ny = 5   
name = ""

if len(sys.argv) > 1:
  if "check_results.py" not in sys.argv[0]:
    name = sys.argv[0]
    n = int(sys.argv[1])
    nx = int(np.sqrt(n))
    ny = nx+1
    
    if "fixed" in name:
      ny = nx
    
    print("name: \"{}\", nx,ny: {},{}".format(name, nx, ny))

# boundary conditions
bc = {}
n_nodes_x = 6
n_nodes_y = 6
if "linear" in name or "hermite" in name:
  n_nodes_x = nx + 1
  n_nodes_y = ny + 1

if "quadratic" in name:
  n_nodes_x = 2*nx + 1
  n_nodes_y = 2*ny + 1

stride = 1
if "hermite" in name:
  stride = 4

k = 1  # mode number

for ix in range(int(n_nodes_x)):
  # x position
  x = float(ix)/(n_nodes_x-1)
  
  # bottom
  i = ix
  bc[stride*i] = 0.0
  
  # top
  i = n_nodes_x*(n_nodes_y-1) + ix
  bc[stride*i] = np.sin(k*np.pi*x)
  
  # left
  i = ix*n_nodes_x
  bc[stride*i] = 0.0
  
  # right
  i = ix*n_nodes_x+(n_nodes_x-1)
  bc[stride*i] = 0.0

physicalExtent = [1.0, 1.0]

#if "fixed" in name:
#  physicalExtent = [nx, ny]

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
      
config = {
  "FiniteElementMethod" : {
    "nElements": [nx, ny],
    "nodePositions": node_positions,
    "elements": elements,
    "physicalExtent": physicalExtent,
    "prefactor": 0.6,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "maxIterations": 500000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+name, "binary": "false", "fixedFormat": False},
      {"format": "ExFile", "filename": "out/"+name, "outputInterval": 2},
      {"format": "PythonFile", "filename": "out/"+name, "outputInterval": 2, "binary":False, "onlyNodalValues":True},
    ]
  },
}
