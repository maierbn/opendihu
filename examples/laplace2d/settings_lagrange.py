import numpy as np

n = 4.
m = n

# boundary conditions
bc = {}
for i in range(int(2*n+1)):
  x = i/(2*n+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (2*n+1)*(2*m) + i
  bc[i2] = np.sin(x*np.pi)
  
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
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}
