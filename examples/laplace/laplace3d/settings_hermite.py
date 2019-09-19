
nx = 3   # number of elements in x direction
ny = 2   # number of elements in y direction
nz = 4   # number of elements in z direction

# boundary conditions
bc = {}
for i in range(int(nx+1)):
  for j in range(int(ny+1)):
    x = i/(nx+1.)
    y = j/(ny+1.)
    bc[8*(j*(nx+1)+i)] = i

    i2 = nz*(ny+1)*(nx+1) + j*(nx+1)+i
    bc[8*i2] = 10.*i

config = {
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,
    "nElements": [nx, ny, nz],
    "physicalExtent": [2*nx, 2*ny, 2*nz],
    "dirichletBoundaryConditions": bc,
    "maxIterations": 1e5,
    "relativeTolerance": 1e-15,
    "OutputWriter" : [
      {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary": False}
    ]
  }
}
