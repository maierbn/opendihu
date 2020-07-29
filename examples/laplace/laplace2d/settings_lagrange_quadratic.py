# 2D Laplace, quadratic elements

import numpy as np

# number of elements
nx = 10
ny = nx

# number of nodes (for quadratic elements)
mx = 2*nx+1  
my = 2*ny+1

# boundary conditions 
bc = {}
for i in range(mx):
  x = i/mx
  
  # bottom line
  bc[i] = np.sin(x*np.pi)
  
  # top line
  i2 = (my-1)*mx + i
  bc[i2] = np.sin(x*np.pi)
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "laplace",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod": {
    # mesh parameters
    "nElements": [nx, ny],
    "inputMeshIsGlobal": True,
    "physicalExtent": [1.0, 1.0],
    "outputInterval": 1.0,
    
    # problem parameters
    "dirichletBoundaryConditions": bc,
    "neumannBoundaryConditions": [],
    "prefactor": 1,
    
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
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":False, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ]
  },
}
