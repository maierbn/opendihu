# parallel fiber estimation, Laplace 3D
# refine given fiber file and interpolate fibers in between
#
# command arguments: <nFineGridFibers> <input_filename> <bottom_z_clip> <top_z_clip> <element_length>

import numpy as np
import sys
import pickle

print("settings_refine.py args:",sys.argv)

nFineGridFibers = 1
if len(sys.argv) > 0+2:
  nFineGridFibers = (int)(sys.argv[0])
    
input_filename = "7x7fibers.bin"
if len(sys.argv) > 1+2:
  input_filename = sys.argv[1]

bottom_z_clip = 72.0
if len(sys.argv) > 2+2:
  bottom_z_clip = (float)(sys.argv[2])
    
top_z_clip = 220.0
if len(sys.argv) > 3+2:
  top_z_clip = (float)(sys.argv[3])

element_length = 0.1 # [cm]
if len(sys.argv) > 4+2:
  element_length = (float)(sys.argv[4])

if nFineGridFibers == 0:
  print("Error, nFineGridFibers is 0")
  exit

n_nodes_per_fiber = (top_z_clip-bottom_z_clip) / element_length
n_nodes_per_fiber = 2*(n_nodes_per_fiber//2)+1   # make number odd

print("input_filename: \"{}\"".format(input_filename))
print("nFineGridFibers: {}".format(nFineGridFibers))
print("bottom z clip: {}, top z clip: {}".format(bottom_z_clip, top_z_clip))
print("element_length: {}".format(element_length))

bc = {}

config = {
  "logFormat": "csv",
  "solverStructureDiagramFile": "solver_structure.txt",
  "mappingsBetweenMeshesLogFile": "",
  "scenarioName": "",
  "Solvers": {
    "linearSolver": {
      "relativeTolerance": 1e-12,
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations": 500000,
    }
  },
  "ParallelFiberEstimation" : {
    "inputMeshFilename":        "",                 # not relevant here
    "stlFilename":              "",                 # not relevant here
    "resultFilename":           input_filename,     # the filename to work on
    "waitIfFileGetsBig":        False,              # do not wait for key press if the output files is greater than 1 GB
    "bottomZClip":              bottom_z_clip,      # 82 (72), bottom z value of the muscle volume  
    "topZClip":                 top_z_clip,         # 250 (220), top z value of the muscle volume
    "finalBottomZClip":         bottom_z_clip,      # 82 (72), bottom z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "finalTopZClip":            top_z_clip,         # 250 (220), top z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "useNeumannBoundaryConditions": True,           # which type of boundary conditions at top and bottom should be used, Neumann or Dirichlet type  
    "nElementsXPerSubdomain":   4,                  # number of elements in x and y-direction per subdomain
    "nElementsZPerSubdomain":   50,                 # number of elements in z-direction per subdomain
    "nFineGridFibers":          nFineGridFibers,    # number of additional fine fibers that are interpolated between the main "key" fibers, the key fibers are traced
    "useGradientField":         False,              # set to False
    "maxLevel":                 2,                  # maximum level (1=8 processes, 2=64 processes)
    "lineStepWidth":            0.1,                # line width for tracing of fibers
    "nNodesPerFiber":           n_nodes_per_fiber,  # number of nodes in each final fiber
    "improveMesh":              True,               # smooth the 2D meshes, required for bigger meshes or larger amount of ranks
    "refinementFactors":        [1,1,1],            # no refinement
    "maxIterations":            1e5,
    "laplacianSmoothingNIterations": 10,            # number of Laplacian smoothing iterations on the final fibers grid
    
    "FiniteElementMethod" : {
      "meshName":       "potentialFlow",
      "solverName":     "linearSolver",
      "dirichletBoundaryConditions": bc,
      "dirichletOutputFilename":     None,    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "prefactor":      1.0,
    },
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/bin", "binary": True, "fixedFormat": False, "combineFiles": False, "fileNumbering": "incremental"},
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
