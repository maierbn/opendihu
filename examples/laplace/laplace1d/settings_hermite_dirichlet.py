
# Laplace 1D Hermite

n = 5

# Dirichlet boundary conditions
bc = {}
bc[0] = 1.0       # Dirichlet BC for left dof
bc[1] = +1./5     # Dirichlet BC for left dof, derivative, -1/5 would give the straight line
bc[2*n] = 0.0     # Dirichlet BC for right dof
bc[2*n+1] = -1./5 # Dirichlet BC for right dof, derivative

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "laplace",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
    
    # mesh parameters
    "nElements": n,
    "physicalExtent": 5.0,
    "nodePositions": [[float(i)/n] for i in range(n+1)],
    "elements":      [[i, i+1] for i in range(n)],
    "inputMeshIsGlobal": True,
    
    # problem parameters
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    "prefactor": 1.0,
    
    # solver parameters
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 1e4,
    "dumpFormat": "matlab",
    "dumpFilename": "out/hermite",
    "slotName": "solution",                                     # name of the data slot, this is only relevant if multiple solvers are used that share some data
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "combineFiles": False, "onlyNodalValues": True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/p", "combineFiles": False, "binary": True, "onlyNodalValues": False, "fileNumbering": "incremental"}
    ]
  },
}
