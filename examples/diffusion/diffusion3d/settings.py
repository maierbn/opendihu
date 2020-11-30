# Diffusion 3D
# number of elements
nx = 20   
ny = 8
nz = 10

# number of nodes
mx = nx + 1
my = ny + 1
mz = nz + 1

# initial values
iv = {}

for k in range(int(0.3*nz), int(0.5*nz)):
  for j in range(int(0.2*ny), int(0.3*ny)):
    for i in range(int(0.5*nx), int(0.8*nx)):
      index = k*mx*my+ j*mx + i
      iv[index] = 1.0

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "diffusion",                # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "Heun" : {
    "initialValues": iv,
    "timeStepWidth": 1e-2,
    "endTime": 10.0,
    "timeStepOutputInterval": 10,
    "dirichletBoundaryConditions": {},
    "dirichletOutputFilename": None,   # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "inputMeshIsGlobal": True,
    "checkForNanInf": False,
    "nAdditionalFieldVariables": 0,
    "additionalSlotNames": [],
    
    "FiniteElementMethod" : {
      # mesh parameters
      "nElements": [nx,ny,nz],
      "physicalExtent": [4.0*nx,4.0*ny,4.0*nz],
      "physicalOffset": [0, 0, 0],
      "prefactor": 3,
      "inputMeshIsGlobal": True,
      
      # solver parameters
      "solverType": "cg",
      "preconditionerType": "none",
      "relativeTolerance": 1e-5,         # relative tolerance of the residual normal, respective to the initial residual norm, linear solver
      "absoluteTolerance": 1e-5,         # 1e-10 absolute tolerance of the residual    
      "maxIterations": 1000,
      "dumpFormat": "default",
      "dumpFilename": "",
      "slotName": "",
    },
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 10, "filename": "out/out", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"},
      #{"format": "PythonFile", "outputInterval": 10, "filename": "out/out_diffusion2d", "binary": True, "fileNumbering": "incremental"}
    ]
  },
}
