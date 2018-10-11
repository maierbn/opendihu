import numpy as np

n = 20.
m = n

# boundary conditions
bc = {}
for i in range(int(n+1)):
  x = i/(n+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (n+1)*m + i
  bc[(n+1)*m+i] = np.sin(x*np.pi)
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))

config = {
  "FiniteElementMethod" : {
    "nElements": [n, m],
    "physicalExtent": [1.0, 1.0],
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/p"}
    ]
  },
}
