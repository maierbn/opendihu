import numpy as np

# Poisson 2D

n = 5.
m = n

# boundary conditions
bc = {}
for i in range(int(n+1)):
  x = i/(n+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (n+1)*m + i
  bc[(n+1)*m+i] = np.sin(x*np.pi)
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))

# right hand side
rhs = np.zeros(int((n+1)*(m+1)))
for y in range(int(m+1)):
  for x in range(int(n+1)):
    i = int(y*(n+1) + x)
    rhs[i] = np.sin(float(x)/(n+1)*np.pi)*np.sin(float(y)/(m+1)*2*np.pi)
    
    print("rhs[{}] = {}".format(i, rhs[i]))

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [n, m],
    "physicalExtent": [1.0, 1.0],
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "rightHandSide" : list(rhs),
    "OutputWriter" : [
      {"format": "Paraview", "interval": 1, "filename": "out/p", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/p", "binary":True}
    ]
  },
}
