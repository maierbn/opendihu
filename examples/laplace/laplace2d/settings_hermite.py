import numpy as np

# number of elements
nx = 10
ny = nx

# number of nodes (for Hermite elements)
mx = nx+1  
my = ny+1

# boundary conditions
bc = {}
for i in range(mx):
  # top line
  x = i/mx
  bc[4*i] = np.sin(x*np.pi)
  
  # bottom line
  i2 = (my-1)*mx + i
  bc[4*i2] = np.sin(x*np.pi)
  
config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "laplace",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
    # mesh parameters
    "nElements": [nx, ny],
    "physicalExtent": [1.0, 1.0],
    "physicalOffset": [0, 0],
    "outputInterval": 1.0,
    "setHermiteDerivatives": True,
    "inputMeshIsGlobal": True,
    
    # problem parameters
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    "prefactor": 1,
    
    # solver parameters
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 1e6,
    "dumpFilename": "out/",
    "dumpFormat": "matlab",  # default, ascii, or matlab
    "slotName": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"}
    ]
  },
}
