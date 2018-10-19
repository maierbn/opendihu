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

with_aponeurosis = False

# read in data
with open('../mesh', 'rb') as f:
  data = pickle.load(f, encoding='latin1')

# data contains the following entries:
#  "node_positions": node_positions, 
#  "linear_elements": linear_elements, 
#  "quadratic_elements": quadratic_elements, 
#  "seed_points": seed_points,
#  "bottom_nodes": bottom_node_indices,
#  "top_nodes": top_node_indices,
#  "n_linear_elements_per_coordinate_direction": n_linear_elements_per_coordinate_direction,
#  "n_quadratic_elements_per_coordinate_direction": n_quadratic_elements_per_coordinate_direction,

  node_positions = data["node_positions"]
  seed_points = data["seed_points"]
   
# elements
if "quadratic" in name:
  # quadratic elements
  elements = data["quadratic_elements"]
  n_elements = data["n_quadratic_elements_per_coordinate_direction"]
  
else:
  # linear elements
  elements = data["linear_elements"]
  n_elements = data["n_linear_elements_per_coordinate_direction"]

# boundary conditions
bc = {}
if "hermite" in name:
  for bottom_node_index in data["bottom_nodes"]:
    bc[8*bottom_node_index] = 0.0
  for top_node_index in data["top_nodes"]:
    bc[8*top_node_index] = 1.0
else:
  for bottom_node_index in data["bottom_nodes"]:
    bc[bottom_node_index] = 0.0
  
  if with_aponeurosis:
    if "quadratic" in name:
      n_dofs = [(int)(2*n_elements[0]+1), (int)(2*n_elements[1]+1), (int)(2*n_elements[2]+1)]
    else:
      n_dofs = [n_elements[0]+1, n_elements[1]+1, n_elements[2]+1]
    
    # set top BC along aponeurosis
    z_begin = (int)(0.5*n_dofs[2])
    z_end = n_dofs[2]
    for z in range(z_begin, z_end):
      i = z*n_dofs[0]*n_dofs[1] + (int)(n_dofs[1]/2)*n_dofs[0] + (int)(n_dofs[0]/2)
      
      alpha = (z-z_begin) / (z_end-z_begin)
      bc[i] = 1.0 + alpha*0.1 
      
  else:
    for top_node_index in data["top_nodes"]:
      bc[top_node_index] = 1.0
  
# use the gradient field for linear functions
use_gradient_field = ("linear" in name)

config = {
  "Meshes": {
    "potentialFlow": {
      "nodePositions": node_positions,
      "elements": elements,
      "nElements": n_elements
    }
  },
  "Solvers": {
    "linearSolver": {
      "relativeTolerance": 1e-15,
      "maxIterations": 500000,
    }
  },
  "StreamlineTracer" : {
    "seedPoints": seed_points,
    "maxIterations": 500000,
    "useGradientField": use_gradient_field,
    "lineStepWidth": 1e-1,
    "targetElementLength": 1./100.,  # length per element, i.e. distance between nodes, length = 1/100 cm (100 per cm)
    "targetLength": 15,     # length of longest streamline, all streamlines are equally scaled to fit this value
    "discardRelativeLength": 0.7,
    "csvFilename": "{}.csv".format(name),
    "csvFilenameBeforePostprocessing": "{}_raw.csv".format(name),
    "FiniteElementMethod" : {
      "meshName": "potentialFlow",
      "solverName": "linearSolver",
      "dirichletBoundaryConditions": bc,
      "prefactor": 1.0,
    },
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+name, "binary": True, "fixedFormat": False},
      {"format": "ExFile", "filename": "out/"+name, "outputInterval": 2},
      #{"format": "PythonFile", "filename": "out/"+name, "binary":False, "onlyNodalValues":True},
    ]
  }
}

# output config in a readable format
if False:
  import pprint 
  pp = pprint.PrettyPrinter()
  pp.pprint(config)
