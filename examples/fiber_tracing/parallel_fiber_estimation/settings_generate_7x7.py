# parallel fiber estimation, Laplace 3D
#
# command arguments: <name>

import numpy as np
import sys
import pickle

name = ""

if len(sys.argv) < 7:
  print("arguments: <input_filename> <output_filename> <bottom z> <top z> <element size>")
  sys.exit(0)
else:
  input_filename = sys.argv[0]        # e.g. ../../electrophysiology/input/biceps_splines.stl
  output_filename = sys.argv[1]       # e.g. result_0x0fibers.bin
  bottom_z = float(sys.argv[2])
  top_z = float(sys.argv[3])
  element_size = float(sys.argv[4])

n_nodes_per_fiber = (top_z-bottom_z) / element_size
n_nodes_per_fiber = 2*(n_nodes_per_fiber//2)+1   # make number odd

print("input_filename:  {}".format(input_filename))
print("output_filename: {}".format(output_filename))
print("z range to extract: [{}, {}]".format(bottom_z, top_z))
print("element_size: {}, number of elements per fiber: {}".format(element_size, n_nodes_per_fiber))

config = {
  "scenarioName": "7x7fibers",
  "logFormat": "csv",
  "solverStructureDiagramFile": "solver_structure.txt",
  "mappingsBetweenMeshesLogFile": "",
  "Solvers": {
    "linearSolver": {
      "relativeTolerance": 1e-4,
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations": 1e3,
      "solverType": "gmres",
      "preconditionerType": "sor",
      "dumpFilename": "",
      "dumpFormat": "default",
    }
  },
  "ParallelFiberEstimation" : {
    "inputMeshFilename": input_filename,  # this is the input filename
    "resultFilename":    output_filename, # this is the output filename, the numbers <a>x<b> are adjusted automatically
    "bottomZClip":       bottom_z,        # 82 (72), bottom z value of the muscle volume to simulate the potential flow in
    "topZClip":          top_z,           # 250 (220), top z value of the muscle volume
    "finalBottomZClip":  bottom_z,        # 82 (72), bottom z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "finalTopZClip":     top_z,           # 250 (220), top z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "useNeumannBoundaryConditions": True, # which type of boundary conditions at top and bottom should be used, Neumann or Dirichlet type  
    "nElementsXPerSubdomain": 4,          # 4 number of elements in x and y-direction per subdomain
    "nElementsZPerSubdomain": 10,         # number of elements in z-direction per subdomain
    #"nElementsZPerSubdomain": 50,         # number of elements in z-direction per subdomain
    "nFineGridFibers":   0,               # number of additional fine fibers that are interpolated between the main "key" fibers, the key fibers are traced
    "useGradientField":  False,           # set to False
    "maxLevel":          0,               # maximum level (0=1 process, 1=8 processes, 2=64 processes)
    "lineStepWidth":     0.01,            # line width for tracing of fibers
    "nNodesPerFiber":    n_nodes_per_fiber,   # number of nodes in each final fiber
    "improveMesh":       True,            # smooth the 2D meshes, required for bigger meshes or larger amount of ranks
    #"refinementFactors": [2,2,2],        # [2,2,2] factors in x,y,z direction by which the mesh should be refined prior to solving the laplace problem and tracing the streamlines
    "refinementFactors": [1,1,1],         # no refinement
    "maxIterations":     1e5,             # maximum number of iterations per fiber
    "FiniteElementMethod" : {
      "meshName": "potentialFlow",
      "solverName": "linearSolver",
      "dirichletBoundaryConditions": {},
      "dirichletOutputFilename": None,    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "prefactor": 1.0,
      "inputMeshIsGlobal": True,
      "slotName": "",
    },
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/biceps", "binary": True, "fixedFormat": False, "combineFiles": False, "fileNumbering": "incremental"},
      #{"format": "Paraview", "outputInterval": 1, "filename": "out/txt", "binary": False, "fixedFormat": False, "combineFiles": False},
      #{"format": "ExFile", "filename": "out/"+name, "outputInterval": 2, "sphereSize": "0.005*0.005*0.01"},
      #{"format": "PythonFile", "filename": "out/"+name, "binary":False, "onlyNodalValues":True},
    ]
  }
}
