import numpy as np

# Poisson 2D

n = 5.
m = n

# boundary conditions
bc = {}
for i in range(int(n+1)):
  x = i/(n+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (n+1)*m + i
  bc[(n+1)*m+i] = np.sin(x*np.pi)
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))

# right hand side
rhs = np.zeros(int((n+1)*(m+1)))
for y in range(int(m+1)):
  for x in range(int(n+1)):
    i = int(y*(n+1) + x)
    rhs[i] = np.sin(float(x)/(n+1)*np.pi)*np.sin(float(y)/(m+1)*2*np.pi)
    
    print("rhs[{}] = {}".format(i, rhs[i]))

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "poisson",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
    # mesh parameters
    "nElements": [n, m],
    "physicalExtent": [1.0, 1.0],
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
    "maxIterations": 1000,
    "dumpFormat": "ascii",
    "dumpFilename": "out/poisson",
    "slotName": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": "false", "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
    ]
  },
}
