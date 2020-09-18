# parallel fiber estimation, Laplace 3D
#
# command arguments: splines_or_stl refinement improve_mesh use_gradient_field use_neumann_bc

import numpy as np
import sys
import pickle

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

if len(sys.argv) >= 5:
  splines_or_stl = sys.argv[0]
  refinement = (int)(sys.argv[1])
  improve_mesh = True if sys.argv[2] == "true" else False
  use_gradient_field = True if sys.argv[3] == "true" else False
  use_neumann_bc = True if sys.argv[4] == "true" else False

  scenario_name = "72-220 {} {} {} {}".format(splines_or_stl,refinement,"I" if improve_mesh else "ni","G" if use_gradient_field else "ng","neumann" if use_neumann_bc else "dirichlet")
  result_filename = "0x0fibers_{}{}{}{}{}_out.bin".format(splines_or_stl,refinement,"I" if improve_mesh else "n","G" if use_gradient_field else "n", "_neumann" if use_neumann_bc else "_dirichlet")
  if splines_or_stl == "splines":
    input_mesh_filename = "../../../electrophysiology/input/biceps.surface.pickle"
  else:
    input_mesh_filename = "../../../electrophysiology/input/biceps_splines.stl"

  if rank_no == 0:
    print("parameters:")
    print("splines_or_stl:      {}".format(splines_or_stl))
    print("refinement:          {}".format(refinement))
    print("improve_mesh:        {}".format(improve_mesh))
    print("use_gradient_field:  {}".format(use_gradient_field))
    print("use_neumann_bc:      {}".format(use_neumann_bc))
    print("input_mesh_filename: {}".format(input_mesh_filename))
    print("result_filename:     {}".format(result_filename))

else:
  if rank_no == 0:
    print("error, usage: splines_or_stl refinement improve_mesh use_gradient_field use_neumann_bc")
  sys.exit(0)
  
bc = {}

config = {
  "scenarioName": scenario_name,
  "solverStructureDiagramFile": None,
  "mappingsBetweenMeshesLogFile": None,
  "logFormat": "csv",
  "Solvers": {
    "linearSolver": {
      "relativeTolerance": 1e-4,
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations": 1e3,
      "solverType": "gmres",
      "preconditionerType": "sor",
      "dumpFormat": None,
      "dumpFilename": None,
    }
  },
  "ParallelFiberEstimation" : {
    #"inputMeshFilename": "../../../../testing/system_testing/tests/fibers/meshes/biceps_full.stl",   # this is the input filename
    #"inputMeshFilename": "../../../electrophysiology/input/biceps_full.stl",   # this is the input filename
    #"inputMeshFilename": "../../../electrophysiology/input/biceps_splines.stl",   # this is the input filename
    #"inputMeshFilename": "../../../electrophysiology/input/biceps.surface.pickle",   # this is the input filename
    "inputMeshFilename":          input_mesh_filename,   # this is the input filename
    "resultFilename":             result_filename,       # this is the output filename, the numbers <a>x<b> are adjusted automatically
    
    "bottomZClip":                72.0,         # 82 (72), bottom z value of the muscle volume to simulate the potential flow in
    "topZClip":                   220.0,        # 250 (220), top z value of the muscle volume
    "finalBottomZClip":           72.0,         # 82 (72), bottom z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "finalTopZClip":              220.0,        # 250 (220), top z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "useNeumannBoundaryConditions": use_neumann_bc,      # which type of boundary conditions at top and bottom should be used, Neumann or Dirichlet type  
    "nElementsXPerSubdomain":     4,            # 4 number of elements in x and y-direction per subdomain
    "nElementsZPerSubdomain":     50,           # number of elements in z-direction per subdomain
    "nFineGridFibers":            0,            # number of additional fine fibers that are interpolated between the main "key" fibers, the key fibers are traced
    "useGradientField":           use_gradient_field,    # set to False
    "maxLevel":                   2,            # maximum level (0=1 process, 1=8 processes, 2=64 processes)
    "lineStepWidth":              0.01,         # line width for tracing of fibers
    "nNodesPerFiber": (220.-72.) / 0.1,         # number of nodes in each final fiber
    "maxIterations":              1e5,
    
    "improveMesh": improve_mesh,                # smooth the 2D meshes, required for bigger meshes or larger amount of ranks
    "refinementFactors": [refinement,refinement,refinement],         # [2,2,2] factors in x,y,z direction by which the mesh should be refined prior to solving the laplace problem and tracing the streamlines
    
    "FiniteElementMethod": {
      "meshName": "potentialFlow",
      "solverName": "linearSolver",
      "dirichletBoundaryConditions": bc,
      "dirichletOutputFilename":  None,    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "prefactor": 1.0,
      "maxIterations": 1e5,
      "inputMeshIsGlobal": True,
      "neumannBoundaryConditions": [],
      "slotName": "",
    },
    "OutputWriter": [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/biceps", "binary": True, "fixedFormat": False, "combineFiles": False, "fileNumbering": "incremental"},
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

#sys.exit(0)
