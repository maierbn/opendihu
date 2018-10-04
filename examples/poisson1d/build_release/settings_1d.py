import numpy as np

# Poisson 1D

n = 50


# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

# right hand side
rhs = np.zeros(n+1)
for i in range(n+1):
  x = float(i)/(n+1)
  rhs[i] = np.sin(x*np.pi)
  if i < (n+1)/2:
    rhs[i] = -1
  else:
    rhs[i] = 0
  print("  rhs[{}] = {}".format(i, rhs[i]))

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "FiniteElementMethod" : {
    "nElements": [n],
    "physicalExtent": [1.0],
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "rightHandSide" : list(rhs),
    
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out/p", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/p"}
    ]
  },
}
