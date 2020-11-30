# Laplace 2D, Neumann BC


# number of nodes
nx = 6
ny = nx

# set neumann boundary conditions
neumann_bc = []
import numpy as np

# bottom, top
if False:
  for i in range(nx):
    neumann_bc += [
        {"element": i, "constantValue": -1.0, "face": "1-"},               # bottom
        {"element": (ny-1)*nx + i, "constantValue": 1.0, "face": "1+"},   # top
      ]

# left, right
if True:
  for j in range(ny):
    neumann_bc += [
        {"element": j*nx, "constantValue": -np.sin(j/ny*np.pi), "face": "0-"},            # left
        {"element": j*nx + (nx-1), "constantValue": np.sin(j/ny*np.pi), "face": "0+"},   # top
      ]
      
config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "laplace",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
    "nElements": [nx, ny],
    "inputMeshIsGlobal": True,
    "physicalExtent": [1.0, 1.0],
    "nodePositions": [[float(i)/nx, float(j)/ny, 0] for j in range(ny+1) for i in range(nx+1)],
    "setHermiteDerivatives": True,
    "elements": [[j*(nx+1)+i, j*(nx+1)+i+1, (j+1)*(nx+1)+i, (j+1)*(nx+1)+i+1] for j in range(ny) for i in range(nx)],
    
    #"dirichletBoundaryConditions": {
    #  0:0, 1:0.2, 2:0.4, 3:0.6, 4:0.8, 5:1,   
    #  30:0, 31:0.2, 32:0.4, 33:0.6, 34:0.8, 35:1
    #},
    "dirichletBoundaryConditions": {0:0},
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": neumann_bc,        # Neumann BC are always interpolated using Lagrange ansatz functions with one dof per node (not Hermite), even if the solution uses Hermite ansatz functions
    "prefactor": 1,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 10000,
    "dumpFilename": "out/",
    "dumpFormat": "matlab",  # default, ascii, or matlab
    "slotName": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":False, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ]
  }
}

print(config)
print(len(config["FiniteElementMethod"]["nodePositions"]))
