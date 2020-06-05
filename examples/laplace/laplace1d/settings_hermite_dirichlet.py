
# Laplace 1D Hermite

n = 5

# Dirichlet boundary conditions
bc = {}
bc[0] = 1.0
bc[1] = -1./5
bc[2*n] = 0.0
bc[2*n+1] = -1./5

config = {
  "logFormat":  "csv",
  "solverStructureDiagramFile":     "solver_structure.txt",
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 5.0,
    "nodePositions": [[float(i)/n] for i in range(n+1)],
    "elements": [[i, i+1] for i in range(n)],
    "inputMeshIsGlobal": True,
    
    "dirichletBoundaryConditions": bc,
    "neumannBoundaryConditions": [],
    "prefactor": 1.0,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 1e4,
    "dumpFormat": "matlab",
    "dumpFilename": "out/hermite",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "combineFiles": False, "onlyNodalValues": True},
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/p", "combineFiles": False, "binary": True, "onlyNodalValues": False}
    ]
  },
}
