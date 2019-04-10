
# Laplace 1D

n = 40


# Dirichlet boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 4.0,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "inputMeshIsGlobal": True,
    "prefactor": 1.0,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "combineFiles": False},
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/p", "combineFiles": False}
    ]
  },
}
