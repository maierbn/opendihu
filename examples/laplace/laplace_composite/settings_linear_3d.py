import numpy as np
import sys

n_ranks = (int)(sys.argv[-1])

nx1 = 1
nx2 = 1
ny1 = 1
ny2 = 1
nz1 = 3
nz2 = 2

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "Meshes": {
    "submesh0": {
      "nElements": [nx1, ny1, nz1],
      "inputMeshIsGlobal": True,
      "physicalExtent": [nx1, ny1, nz1],
    },
    "submesh1": {
      "nElements": [nx2, ny2, nz2],
      "inputMeshIsGlobal": True,
      "physicalExtent": [nx2, ny2, nz2],
      "physicalOffset": [0.0, ny1, 0.0],
    },
  },
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,
    "outputInterval": 1.0,
    
    "dirichletBoundaryConditions": {0:1,-1:-1},
    "neumannBoundaryConditions": [],
    "prefactor": [1,2],
    
    "meshName": ["submesh0", "submesh1"],
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "maxIterations": 10000,
    "dumpFormat": "default",
    "dumpFilename": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/linear_3d_{}".format(n_ranks), "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      #{"format": "PythonFile", "filename": "out/2d", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}
