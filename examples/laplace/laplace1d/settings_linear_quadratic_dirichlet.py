
# Laplace 1D

n = 40


# Dirichlet boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "FiniteElementMethod" : {
  
    "nElements": n,
    "inputMeshIsGlobal": True,
    "physicalExtent": 4.0,
  
    "dirichletBoundaryConditions": bc,
    "neumannBoundaryConditions": [],
    "prefactor": 1.0,
  
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "maxIterations": 1e4,
    "dumpFormat": "matlab",
    "dumpFilename": "out/dirichlet",
  
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "combineFiles": False, "onlyNodalValues": True},
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/p", "binary": False, "combineFiles": False, "onlyNodalValues": True}
    ]
  },
}
