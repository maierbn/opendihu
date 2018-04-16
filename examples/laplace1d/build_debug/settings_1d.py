
# Laplace 1D

n = 40


# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 4.0,
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/p"}
    ]
  },
}
