# parallel fiber estimation, Laplace 3D
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

bc = {}

config = {
  "Solvers": {
    "linearSolver": {
      "relativeTolerance": 1e-15,
      "maxIterations": 500000,
    }
  },
  "ParallelFiberEstimation" : {
    "stlFilename": "../../../testing/system_testing/tests/fibers/meshes/biceps_full.stl",
    "bottomZClip":  37.0,   # top z value of the muscle volume
    "topZClip": 300.0,      # bottom z value of the muscle volume
    "nElementsZPerSubdomain": 4,  # number of elements in z-direction per subdomain
    "maxLevel": 1,          # maximum level
    "FiniteElementMethod" : {
      "meshName": "potentialFlow",
      "solverName": "linearSolver",
      "dirichletBoundaryConditions": bc,
      "prefactor": 1.0,
    },
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+name, "binary": True, "fixedFormat": False, "combineFiles": False},
      {"format": "ExFile", "filename": "out/"+name, "outputInterval": 2, "sphereSize": "0.005*0.005*0.01"},
      #{"format": "PythonFile", "filename": "out/"+name, "binary":False, "onlyNodalValues":True},
    ]
  }
}

# output config in a readable format
if False:
  import pprint 
  pp = pprint.PrettyPrinter()
  pp.pprint(config)
