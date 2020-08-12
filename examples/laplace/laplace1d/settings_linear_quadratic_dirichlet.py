
# Laplace 1D

n = 40


# Dirichlet boundary conditions
bc = {}
bc[0] = 1.0         # set degree of freedom 0 (left most) to 1
bc[-1] = 0.0        # -1 means the last degree of freedom, set right most degree of freedom to 0

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "laplace",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
  
    # mesh parameters
    "nElements":          n,
    "inputMeshIsGlobal":  True,
    "physicalExtent":     4.0,
    "physicalOffset":     0,
  
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
    "maxIterations":      1e4,
    "dumpFormat": "matlab",
    "dumpFilename": "out/dirichlet",
    "slotName": "solution",                                     # name of the data slot, this is only relevant if multiple solvers are used that share some data
  
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "combineFiles": False, "onlyNodalValues": True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/p", "binary": False, "combineFiles": False, "onlyNodalValues": True, "fileNumbering": "incremental"}
    ]
  },
}
