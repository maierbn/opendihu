
# Laplace 1D

n = 40


# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "numCycles" : 100,
  "Term1":{
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 4.0,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "Paraview", "interval": 1, "filename": "out/p", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/p"}
    ]
  }
  },
  "Term2":{
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 4.0,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "Paraview", "interval": 1, "filename": "out/p1", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/p2"}
    ]
  }
  },
}
