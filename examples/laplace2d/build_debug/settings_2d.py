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
    "outputInterval": 1.0,
    "prefactor": 1,
    "DirichletBoundaryCondition": bc,
    "relativeTolerance": 1e-15,
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binaryOutput": "false", "fixedFormat": False, "onlyNodalValues":True},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}
