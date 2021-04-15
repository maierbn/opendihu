import numpy as np

# Poisson 1D
n = 11

# Dirichlet boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

# right hand side
rhs = np.zeros(n+1)
for i in range(n+1):
  x = float(i)/(n+1)
  rhs[i] = np.sin(x*np.pi)
  if i < int(n/2):
    rhs[i] = -10
  else:
    rhs[i] = 0
  print("  rhs[{}] = {}".format(i, rhs[i]))

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "poisson",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
    # mesh parameters
    "nElements": [n],
    "physicalExtent": [1.0],
    "inputMeshIsGlobal": True,
    
    # problem parameters
    "prefactor": 1,
    "rightHandSide": list(rhs),         # provide the rhs vector as list, not as numpy array
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename": None,    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    
    # solver parameters
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
    "maxIterations": 10000,
    "dumpFormat": "ascii",
    "dumpFilename": "out/poisson",
    "slotName": "",
    
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out/p", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
    ]
  },
}
