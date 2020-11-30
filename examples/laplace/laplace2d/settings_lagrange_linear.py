import numpy as np

n = 10
m = n

# boundary conditions (for linear elements)
bc = {}
for i in range(int(n+1)):
  x = i/(n+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (n+1)*m + i
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
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions":   [],
    "prefactor": 1,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "dumpFilename": "out/",
    "dumpFormat": "matlab",  # default, ascii, or matlab
    "maxIterations": 10000,
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":False},
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}
