import numpy as np

# boundary conditions
bc = {}
bc[0] = 0
bc[-1] = 1

config = {
  "FiniteElementMethod" : {
    "nElements": [4,4],
    "physicalExtent": [1.0, 1.0],
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "ExFile", "outputInterval": 1, "filename": "out/p"},
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}
