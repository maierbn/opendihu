import numpy as np

n = 10
m = n

# boundary conditions (for linear elements)
bc = {}
for i in range(int(n+1)):
  x = i/(n+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (n+1)*m + i
  bc[i2] = np.sin(x*np.pi)
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))

config = {
  "FiniteElementMethod" : {
    "nElements": [n, m],
    "inputMeshIsGlobal": True,
    "physicalExtent": [1.0, 1.0],
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "none",
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":False},
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}
