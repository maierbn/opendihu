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
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "laplace_composite",        # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "Meshes": {
    "submesh0": {
      "nElements": [nx1, ny1, nz1],
      "inputMeshIsGlobal": True,
      "physicalExtent": [nx1, ny1, nz1],
      "physicalOffset": [0.0, 0.0, 0.0],
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
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    "prefactor": [1,2],
    
    "meshName": ["submesh0", "submesh1"],
    
    # solver parameters
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 10000,
    "dumpFormat": "default",
    "dumpFilename": "",
    "slotName": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/linear_3d_{}".format(n_ranks), "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      #{"format": "PythonFile", "filename": "out/2d", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ]
  },
}
