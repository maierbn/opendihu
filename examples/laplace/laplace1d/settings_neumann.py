# Laplace 1D
n = 5

# Neumann boundary conditions
bc = [
  {"element": 0, "constantValue": -1.0, "face": "0-"},
  {"element": n-1, "constantValue": 1.0, "face": "0+"}
]
config = {
  "logFormat":                      "csv",
  "solverStructureDiagramFile":     "solver_structure.txt",
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 1.0,
    "nodePositions": [[float(i)/(2*n)] for i in range(2*n+1)],
    "elements": [[2*i, 2*i+1, 2*i+2] for i in range(n)],
    "inputMeshIsGlobal": True,
    
    "neumannBoundaryConditions": bc,
    "dirichletBoundaryConditions": {1:0},
    "prefactor": 1.0,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 1e4,
    "dumpFormat": "matlab",
    "dumpFilename": "out/neumann",
  
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/neumann", "binary": False, "fixedFormat": False, "combineFiles": False, "onlyNodalValues": True},
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/neumann", "combineFiles": False, "binary": True, "onlyNodalValues": True}
    ]
  },
}
