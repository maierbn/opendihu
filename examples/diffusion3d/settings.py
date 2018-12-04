# Diffusion 3D
nx = 40   # number of elements
ny = 20
nz = 30

# initial values
iv = {}

for z in range(int(0.3*nz), int(0.5*nz)):
  for y in range(int(0.2*ny), int(0.3*ny)):
    for x in range(int(0.5*nx), int(0.8*nx)):
      i = z*(nx+1)*(ny+1) + y*(nx+1) + x
      iv[i] = 1.0

#print("iv: ",iv)

config = {
  "Heun" : {
    "initialValues": iv,
    "timeStepWidth": 1e-3,
    "endTime": 1e-1,
    "timeStepOutputInterval": 100,
    "dirichletBoundaryConditions": {},
    "inputMeshIsGlobal": True,
    
    "FiniteElementMethod" : {
      "nElements": [nx,ny,nz],
      "physicalExtent": [4.0*nx,4.0*ny,4.0*nz],
      "relativeTolerance": 1e-15,
      "prefactor": 0.1,
      "inputMeshIsGlobal": True,
      "maxIterations": 10000,
      "solverType": "gmres",
      "preconditionerType": "none",
    },
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/out", "binary": True, "fixedFormat": False, "frequency": 100, "combineFiles": True},
      #{"format": "PythonFile", "outputInterval": 10, "filename": "out/out_diffusion2d", "binary": True}
    ]
  },
}
