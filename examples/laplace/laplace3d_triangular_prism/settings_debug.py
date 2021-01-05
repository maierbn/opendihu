import numpy as np
import sys

# Dirichlet boundary conditions
bc = {}
rank_no = (int)(sys.argv[-2])

# 2
# 0-1

node_positions = [[0,0,0], [1,0,0], [0,1,0], [0.5,0.5,0],
                  [0,0,1], [1,0,1], [0,1,1], [0.5,0.5,1]]

config = {
  "scenarioName": "",
  "mappingsBetweenMeshesLogFile": None,
  "logFormat":                      "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "FiniteElementMethod" : {
    "nElements":  [1,1,1],
    "nRanks":     [1,1,1],
    "hasTriangleCorners": True,
    "inputMeshIsGlobal": True,
    "nodePositions": node_positions,
    "outputInterval": 1,
    "slotName": None,
    
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions":   [],
    "prefactor": 1,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 10000,
    "dumpFormat": "matlab",
    "dumpFilename": "out/a",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},      
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"}
    ]
  },
}
