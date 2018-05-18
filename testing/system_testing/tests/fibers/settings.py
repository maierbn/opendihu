# Laplace 3D
#
# command arguments: <name>

import numpy as np
import sys
import pickle

name = ""

if len(sys.argv) > 0:
  if "check_results.py" not in sys.argv[0]:
    name = sys.argv[0]
    
    print("name: \"{}\"".format(name))

# read in elements and nodes
with open('../mesh', 'rb') as f:
  config_from_file = pickle.load(f, encoding='latin1')

  node_positions = []
  for node_position in config_from_file["nodePositions"]:
    node_positions.append(list(node_position))

# boundary conditions and seed points
bc = {}
seed_points = []

min_z = np.min(node_positions,axis=0)[2]
max_z = np.max(node_positions,axis=0)[2]
print("min z: {}, max z: {}".format(min_z, max_z))

for i,node_position in enumerate(node_positions):
  z = node_position[2]
  
  if z < min_z + 1:
    print("bottom node {}".format(node_position))
    bc[i] = 0.0
    seed_points.append(node_position)
    
  if z > max_z - 1:  
    print("top node {}".format(node_position))
    bc[i] = 1.0



config = {
  "StreamlineTracer" : {
    "seedPoints": seed_points,
    "lineStepWidth": 1e-2,
    "FiniteElementMethod" : {
      "nodePositions": node_positions,
      "elements": config_from_file["elements"],
      "DirichletBoundaryCondition": bc,
      "relativeTolerance": 1e-15,
      "prefactor": 1.0,
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      {"format": "ExFile", "filename": "out/"+name, "outputInterval": 2},
      {"format": "PythonFile", "filename": "out/"+name, "binary":False, "onlyNodalValues":True},
    ]
  }
}

# output config in a readable format
if False:
  import pprint 
  pp = pprint.PrettyPrinter()
  pp.pprint(config)
