# Diffusion 2D
nx = 20   # number of elements
ny = nx

# number of nodes
mx = nx + 1
my = ny + 1

# initial values
iv = {}

for j in range(int(0.2*my), int(0.3*my)):
  for i in range(int(0.5*mx), int(0.8*mx)):
    index = j*mx + i
    iv[index] = 1.0

print("iv: ",iv)

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "diffusion",                # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "Heun" : {
    "initialValues": iv,
    "timeStepWidth": 1e-3,
    "endTime": 1.0,
    "timeStepOutputInterval": 100,
    "checkForNanInf": False,
    "inputMeshIsGlobal": True,
    "nAdditionalFieldVariables": 0,
    "additionalSlotNames": [],
    "dirichletBoundaryConditions": {},
    "dirichletOutputFilename": None,   # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    "FiniteElementMethod" : {
      "nElements": [nx, ny],
      "physicalExtent": [4.0,4.0],
      "inputMeshIsGlobal": True,
      "prefactor": 0.1,
      
      # solver parameters
      "solverType": "gmres",
      "preconditionerType": "none",
      "relativeTolerance": 1e-15,         # relative tolerance of the residual normal, respective to the initial residual norm, linear solver
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations": 10000,
      "dumpFilename": "",
      "dumpFormat": "ascii",  # ascii, default or matlab
      "slotName": "",
    },
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 10, "filename": "out/filename", "binary": "false", "fixedFormat": False, "onlyNodalValues": True, "combineFiles": True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "outputInterval": 10, "filename": "out/out_diffusion2d", "binary": True, "onlyNodalValues": True, "combineFiles": True, "fileNumbering": "incremental"}
    ]
  },
}
