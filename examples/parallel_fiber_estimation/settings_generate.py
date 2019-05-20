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
    
    #print("name: \"{}\"".format(name))

bc = {}

config = {
  "Solvers": {
    "linearSolver": {
      "relativeTolerance": 1e-4,
      "maxIterations": 1e3,
      "solverType": "gmres",
      "preconditionerType": "sor"
    }
  },
  "ParallelFiberEstimation" : {
    "stlFilename": "../../electrophysiology/input/biceps_splines.stl",   # this is the input filename
    "resultFilename": "result_0x0fibers.bin",              # this is the output filename, the numbers <a>x<b> are adjusted automatically
    "bottomZClip":  72.0,                 # 82 (72), bottom z value of the muscle volume to simulate the potential flow in
    "topZClip": 220.0,                    # 250 (220), top z value of the muscle volume
    "finalBottomZClip":  72.0,            # 82 (72), bottom z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "finalTopZClip": 220.0,               # 250 (220), top z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "useNeumannBoundaryConditions": True, # which type of boundary conditions at top and bottom should be used, Neumann or Dirichlet type  
    "nElementsXPerSubdomain": 4,          # number of elements in x and y-direction per subdomain
    "nElementsZPerSubdomain": 50,         # number of elements in z-direction per subdomain
    "nFineGridFibers": 0,                 # number of additional fine fibers that are interpolated between the main "key" fibers, the key fibers are traced
    "useGradientField": False,            # set to False
    "maxLevel": 0,                        # maximum level (0=1 process, 1=8 processes, 2=64 processes)
    "lineStepWidth":  0.1,                # line width for tracing of fibers
    "nNodesPerFiber": (220.-72.) / 0.1,   # number of nodes in each final fiber
    "improveMesh": True,                  # smooth the 2D meshes, required for bigger meshes or larger amount of ranks
    "refinementFactors": [2,2,2],         # factors in x,y,z direction by which the mesh should be refined prior to solving the laplace problem and tracing the streamlines
    "FiniteElementMethod" : {
      "meshName": "potentialFlow",
      "solverName": "linearSolver",
      "dirichletBoundaryConditions": bc,
      "prefactor": 1.0,
    },
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/biceps", "binary": True, "fixedFormat": False, "combineFiles": False},
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
