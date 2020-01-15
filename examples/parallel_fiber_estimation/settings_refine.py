# parallel fiber estimation, Laplace 3D
# refine given fiber file and interpolate fibers in between
#
# command arguments: <nFineGridFibers> <input_filename>

import numpy as np
import sys
import pickle

print("args:",sys.argv)

nFineGridFibers = 1
if len(sys.argv) > 0+2:
  nFineGridFibers = (int)(sys.argv[0])
    
input_filename = "7x7fibers.bin"
if len(sys.argv) > 1+2:
  input_filename = sys.argv[1]
    
if nFineGridFibers == 0:
  print("Error, nFineGridFibers is 0")
  exit

print("input_filename: \"{}\"".format(input_filename))
print("nFineGridFibers: {}".format(nFineGridFibers))

bc = {}

config = {
  "Solvers": {
    "linearSolver": {
      "relativeTolerance": 1e-12,
      "maxIterations": 500000,
    }
  },
  "ParallelFiberEstimation" : {
    "inputMeshFilename": "",  # not relevant here
    "stlFilename": "",        # not relevant here
    "resultFilename": input_filename,
    "waitIfFileGetsBig": False,
    "bottomZClip":  72.0,   # 82 (72), bottom z value of the muscle volume  
    "topZClip": 220.0,      # 250 (220), top z value of the muscle volume
    "finalBottomZClip":  72.0,            # 82 (72), bottom z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "finalTopZClip": 220.0,               # 250 (220), top z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "useNeumannBoundaryConditions": True, # which type of boundary conditions at top and bottom should be used, Neumann or Dirichlet type  
    "nElementsXPerSubdomain": 4,  # number of elements in x and y-direction per subdomain
    "nElementsZPerSubdomain": 50,  # number of elements in z-direction per subdomain
    "nFineGridFibers": nFineGridFibers,     # number of additional fine fibers that are interpolated between the main "key" fibers, the key fibers are traced
    "useGradientField": False,    # set to False
    "maxLevel": 2,          # maximum level (1=8 processes, 2=64 processes)
    "lineStepWidth":  0.1,  # line width for tracing of fibers
    "nNodesPerFiber": (220.-72.) / 0.1,   # number of nodes in each final fiber
    "improveMesh": True,     # smooth the 2D meshes, required for bigger meshes or larger amount of ranks
    "refinementFactors": [1,1,1],         # no refinement
    "maxIterations":     1e5,
    "FiniteElementMethod" : {
      "meshName": "potentialFlow",
      "solverName": "linearSolver",
      "dirichletBoundaryConditions": bc,
      "prefactor": 1.0,
    },
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/bin", "binary": True, "fixedFormat": False, "combineFiles": False},
      #{"format": "Paraview", "outputInterval": 1, "filename": "out/txt", "binary": False, "fixedFormat": False, "combineFiles": False},
      #{"format": "ExFile", "filename": "out/"+name, "outputInterval": 2, "sphereSize": "0.005*0.005*0.01"},
      #{"format": "PythonFile", "filename": "out/"+name, "binary":False, "onlyNodalValues":True},
    ]
  }
}

# output config in a readable format
if False:
  import pprint 
  pp = pprint.PrettyPrinter()
  pp.pprint(config)
