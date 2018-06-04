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

# boundary conditions
bc = {}
bc[0] = 1.0

physicalExtent = 4.0

# set node positions for unstructured mesh
n_nodes = n+1
if "quadratic" in name:
  n_nodes = 2*n + 1
  bc[n_nodes-1] = 0.0
elif "hermite" in name:
  n_dofs = 2*n_nodes
  bc[n_dofs-2] = 0.0
else:
  bc[n] = 0.0  
  
  
node_positions = list(np.linspace(0,4.0,n_nodes))
h = 4.0/(n_nodes)

# move node positions
for i in range(n_nodes):
  node_positions[i] += h*np.sin(float(i)/n_nodes*4*np.pi)

# set elements for unstructured mesh
elements = [[i, i+1] for i in range(n_nodes-1)]

if name == "1d_unstructured_deformable_quadratic":
  elements = [[i, i+1, i+2] for i in range(0,n_nodes-1,2)]

print(node_positions)
print(elements)
print(len(elements))
print(bc)

config = {
  "FiniteElementMethod" : {
    "nElements": n,
    "nodeDimension": 1,
    "nodePositions": node_positions,
    "elements": elements,
    "physicalExtent": physicalExtent,
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+name, "binaryOutput": "false", "fixedFormat": False},
      {"format": "ExFile", "filename": "out/"+name, "outputInterval": 2},
      {"format": "PythonFile", "filename": "out/"+name, "outputInterval": 2, "binary":False, "onlyNodalValues":True},
    ]
  },
}
