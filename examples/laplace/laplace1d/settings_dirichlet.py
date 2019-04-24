
# Laplace 1D Hermite

n = 5

# Dirichlet boundary conditions
bc = {}
bc[0] = 1.0
bc[1] = -1./5
bc[2*n] = 0.0
bc[2*n+1] = -1./5

config = {
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 5.0,
    "dirichletBoundaryConditions": bc,
    "nodePositions": [[float(i)/n] for i in range(n+1)],
    "elements": [[i, i+1] for i in range(n)],
    "relativeTolerance": 1e-15,
    "inputMeshIsGlobal": True,
    "prefactor": 1.0,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "combineFiles": False},
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/p", "combineFiles": False, "onlyNodalValues": False}
    ]
  },
}
