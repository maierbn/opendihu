import numpy as np

n = 4
m = n

# boundary conditions
bc = {}
for i in range(int(n+1)):
  x = i/(n+1.)
  bc[4*i] = np.sin(x*np.pi)
  i2 = int((n+1)*m + i)
  bc[4*i2] = np.sin(x*np.pi)
  
config = {
  "logFormat":                      "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "FiniteElementMethod" : {
    "nElements": [n, m],
    "physicalExtent": [1.0, 1.0],
    "outputInterval": 1.0,
    "setHermiteDerivatives": True,
    "inputMeshIsGlobal": True,
    
    "dirichletBoundaryConditions": bc,
    "neumannBoundaryConditions": [],
    "prefactor": 1,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 1e6,
    "dumpFilename": "out/",
    "dumpFormat": "matlab",  # default, ascii, or matlab
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}
