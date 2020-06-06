import numpy as np

n = 10
m = n

# boundary conditions (for quadratic elements)
bc = {}
for i in range(int(2*n+1)):
  x = i/(2*n+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (2*n+1)*(2*m) + i
  bc[i2] = np.sin(x*np.pi)
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))

config = {
  "logFormat":                      "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "FiniteElementMethod" : {
    "nElements": [n, m],
    "inputMeshIsGlobal": True,
    "physicalExtent": [1.0, 1.0],
    "outputInterval": 1.0,
    
    "dirichletBoundaryConditions": bc,
    "neumannBoundaryConditions": [],
    "prefactor": 1,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
    "maxIterations": 10000,
    "dumpFormat": "default",
    "dumpFilename": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":False},
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}
