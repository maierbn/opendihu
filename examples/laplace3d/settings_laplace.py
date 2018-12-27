import numpy as np

# 3D laplace problem

nx = 10.
ny = nx
nz = nx

# Dirichlet boundary conditions
bc = {}
for j in range(int(ny+1)):
  for i in range(int(nx+1)):
    # x = 0 plane
    x = i/nx
    y = j/ny
    index = int(j*(nx+1) + i)
  
    # z- plane (bottom)
    bc[index] = np.sin(x*np.pi)
    
    # z+ plane (top)
    index += int(nz*(nx+1)*(ny+1))
    bc[index] = np.sin(y*np.pi) + 2.0

config = {
  "FiniteElementMethod" : {
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": [1.0, 1.0, 3.0],
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "none",
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":False},      
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}
