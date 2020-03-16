import numpy as np
import sys

n_ranks = (int)(sys.argv[-1])

nx1 = 4
nx2 = 3
ny = 1

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "Meshes": {
    "submesh0": {
      "nElements": [nx1, ny],
      "inputMeshIsGlobal": True,
      "physicalExtent": [4.0, 1.0],
    },
    "submesh1": {
      "nElements": [nx2, ny],
      "inputMeshIsGlobal": True,
      "physicalExtent": [3.0, 1.0],
      "physicalOffset": [0.0, 1.0],
    },
  },
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,
    "outputInterval": 1.0,
    
    "dirichletBoundaryConditions": {},
    "neumannBoundaryConditions": [],
    "prefactor": [1,2],
    
    "meshName": ["submesh0", "submesh1"],
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 10000,
    "dumpFormat": "default",
    "dumpFilename": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/linear_2d_{}".format(n_ranks), "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      #{"format": "PythonFile", "filename": "out/2d", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}
