import numpy as np

n = 4
m = n

# boundary conditions
bc = {}
for i in range(int(n+1)):
  x = i/(n+1.)
  bc[4*i] = np.sin(x*np.pi)
  i2 = int((n+1)*m + i)
  bc[4*i2] = np.sin(x*np.pi)
  
config = {
  "FiniteElementMethod" : {
    "nElements": [n, m],
    "physicalExtent": [1.0, 1.0],
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "maxIterations": 1e6,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}
