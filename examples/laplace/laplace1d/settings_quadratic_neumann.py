# Laplace 1D
n = 5

# Neumann boundary conditions
bc = [
  {"element": 0, "constantValue": -1.0, "face": "0-"},
  {"element": n-1, "constantValue": 1.0, "face": "0+"}
]
config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "laplace",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   None,                       # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 1.0,
    "nodePositions": [[float(i)/(2*n)] for i in range(2*n+1)],
    "elements": [[2*i, 2*i+1, 2*i+2] for i in range(n)],
    "inputMeshIsGlobal": True,
    
    "neumannBoundaryConditions": bc,
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "dirichletBoundaryConditions": {1:0},
    "prefactor": 1.0,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 1e4,
    "dumpFormat": "matlab",
    "dumpFilename": "out/neumann",
    "slotName": "solution",                                     # name of the data slot, this is only relevant if multiple solvers are used that share some data
  
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/neumann", "binary": False, "fixedFormat": False, "combineFiles": False, "onlyNodalValues": True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/neumann", "combineFiles": False, "binary": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
    ]
  },
}
