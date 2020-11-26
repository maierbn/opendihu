# Laplace 3D
#
# Simulate 3D laplace potential flow and trace streamlines.
# All input data is provided in a pickle file which can be given as <input_filename>
#
# The input pickle file must contain a dict with the following fields:
#  "node_positions": node_positions, 
#  "linear_elements": linear_elements, 
#  "quadratic_elements": quadratic_elements, 
#  "seed_points": seed_points,
#  "bottom_nodes": bottom_node_indices,
#  "top_nodes": top_node_indices,
#  "n_linear_elements_per_coordinate_direction": n_linear_elements_per_coordinate_direction,
#  "n_quadratic_elements_per_coordinate_direction": n_quadratic_elements_per_coordinate_direction,
#
# arguments: <input_filename> <output_filename> [<target_element_length> [<target_fiber_length> [<bottom_to_top>]]]
#
# target_element_length: the mesh width of the fibers in cm
# target_fiber_length:   the length of the longest fiber that is traced, all other fibers are scaled accordingly
# bottom_to_top:         direction, either 0 = top to bottom, or 1 = bottom to top, default: bottom to top

import numpy as np
import sys, os
import pickle
import time

# parse command line arguments
input_filename = "../processed_meshes/mesh_normal"
output_filename = "laplace3d_structured_linear"

target_element_length = 1./100.        # [cm] length per element, i.e. distance between nodes, length = 1/100 cm (100 per cm)
target_fiber_length = 15               # [cm] length of longest streamline, all streamlines are equally scaled to fit this value, set to 0 to disable
bottom_to_top = True

if len(sys.argv)-2 > 0:
  input_filename = sys.argv[0]

if len(sys.argv)-2 > 1:
  output_filename = sys.argv[1]

if len(sys.argv)-2 > 2:
  target_element_length = (float)(sys.argv[2])

if len(sys.argv)-2 > 3:
  target_fiber_length = (float)(sys.argv[3])
    
if len(sys.argv)-2 > 4:
  bottom_to_top = False if (int)(sys.argv[3]) == 0 else True
    
print("input filename: \"{}\"".format(input_filename))
print("output filename: \"{}\"".format(output_filename))

with_aponeurosis = False

# read in data from pickle file
if not os.path.isfile(input_filename):
  time.sleep(1)
with open(input_filename, 'rb') as f:
  data = pickle.load(f, encoding='latin1')

# load node positions and seed points
node_positions = data["node_positions"]
seed_points = data["seed_points"]
 
# create elements
if "quadratic" in output_filename:
  # quadratic elements
  elements = data["quadratic_elements"]
  n_elements = data["n_quadratic_elements_per_coordinate_direction"]
  
else:
  # linear elements
  elements = data["linear_elements"]
  n_elements = data["n_linear_elements_per_coordinate_direction"]

# create boundary conditions
if not bottom_to_top:
  data["bottom_nodes"],data["top_nodes"] = data["top_nodes"],data["bottom_nodes"]
  
bc = {}
if "hermite" in output_filename:
  for bottom_node_index in data["bottom_nodes"]:
    bc[8*bottom_node_index] = 0.0
  for top_node_index in data["top_nodes"]:
    bc[8*top_node_index] = 1.0
else:
  for bottom_node_index in data["bottom_nodes"]:
    bc[bottom_node_index] = 0.0
  
  if with_aponeurosis:
    if "quadratic" in output_filename:
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
use_gradient_field = ("linear" in output_filename)

config = {
  "scenarioName": "streamline_tracer",
  "logFormat": "csv",
  "mappingsBetweenMeshesLogFile": None,
  "solverStructureDiagramFile": None,
  "Meshes": {
    "potentialFlow": {
      "nodePositions": node_positions,
      "elements":      elements,
      "nElements":     n_elements,
      "inputMeshIsGlobal": True,
    }
  },
  "Solvers": {
    "linearSolver": {
      "solverType":         "lu",
      "preconditionerType": "none",
      "relativeTolerance":  1e-15,
      "absoluteTolerance":  0,
      "maxIterations":      5e3,
      "dumpFormat":         "",
      "dumpFilename":       "",
    }
  },
  "StreamlineTracer": {
    "seedPoints":           seed_points,                  # a list of seed points where to start the fibers
    "maxIterations":        5e5,                          # the maximum number of iterations for tracing a single fiber
    "useGradientField":     use_gradient_field,           # if the computed gradient field should be used, There are 2 implementations of streamline tracing. The first one (useGradientField_) uses a precomputed gradient field that is interpolated linearly and the second uses the gradient directly from the Laplace solution field. // The first one seems more stable, because the gradient is zero and the position of the boundary conditions.
    "lineStepWidth":        1e-1,                         # the step width of the tracing
    "targetElementLength":  target_element_length,        # length per element, i.e. distance between nodes, length = 1/100 cm (100 per cm)
    "targetLength":         target_fiber_length,          # length of longest streamline, all streamlines are equally scaled to fit this value, set to 0 to disable
    "discardRelativeLength": 0.7,                         # a relative length (in [0,1]), at the end streamlines are dropped that are smaller than this relative length times the median fiber length
    "csvFilenameBeforePostprocessing": "{}_raw.csv".format(output_filename),   # output file before filtering
    "csvFilename":                     "{}.csv".format(output_filename),       # output file after filtering, contains only valid fibers
    
    # the finite element method object that computes the Laplace equation
    "FiniteElementMethod": {
      "meshName":   "potentialFlow",
      "solverName": "linearSolver",
      "dirichletBoundaryConditions": bc,
      "dirichletOutputFilename":  None,    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "neumannBoundaryConditions": [],
      "prefactor":  1.0,
      "inputMeshIsGlobal": True,
      "slotName": "",
    },
    
    "OutputWriter": [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+output_filename, "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
      #{"format": "ExFile", "filename": "out/"+output_filename, "outputInterval": 2},
      #{"format": "PythonFile", "filename": "out/"+output_filename, "binary":False, "onlyNodalValues":True},
    ]
  }
}

# print the mesh
#print("mesh: {}".format(config["Meshes"]["potentialFlow"]))
