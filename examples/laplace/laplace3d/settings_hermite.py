
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
  "logFormat":                      "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,
    "nElements": [nx, ny, nz],
    "physicalExtent": [2*nx, 2*ny, 2*nz],
    
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    "prefactor": 1,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 10000,
    "dumpFormat": "default",
    "dumpFilename": "",
    
    "OutputWriter" : [
      {"format": "PythonFile", "filename": "out", "outputInterval": 1, "binary": False}
    ]
  }
}
