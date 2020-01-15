import numpy as np

nx = 4
ny = 4

# boundary conditions (for linear elements)
dirichlet_bc = {0: [0.0,0.0,0.0]}

for j in range(1,ny+1):
  dirichlet_bc[j*(nx+1)] = [0.0,None,None]

dirichlet_bc[1] = [0.0,0.0]

neumann_bc = [{"element": j*nx+(nx-1), "constantVector": [0.1,+0.2], "face": "0+"} for j in range(ny)]

#dirichlet_bc = {}
#neumann_bc = []

config = {
  "FiniteElementMethod" : {
    "nElements": [nx, ny],
    "inputMeshIsGlobal": True,
    "physicalExtent": [nx, ny],
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "none",
    "maxIterations": 1e4,
    "bulkModulus": 1.5,
    "shearModulus": 2.0,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}
