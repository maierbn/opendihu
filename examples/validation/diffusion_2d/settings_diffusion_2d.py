import numpy as np

# problem definition
# boundary conditions
def u_boundary(x, y):
  return 1 + x/3 + y/8

# initial values
def u_initial(x, y):
  return np.maximum(0.0, 2 - np.sin(np.sqrt((x-2)**2 + (y-2)**2)))

# discretization
nx = 10   # number of elements
ny = nx

# number of nodes
import opendihu
if "linear" in opendihu.program_name:
  mx = nx + 1
  my = ny + 1

elif "quadratic" in opendihu.program_name:
  mx = 2*nx + 1
  my = 2*ny + 1

# define Dirichlet boundary conditions
bc = {}

# x=0 and x=3
for j in range(my):
  y = 4.0*j/(my-1)
  bc[j*mx + 0] = u_boundary(0,y)
  bc[j*mx + mx-1] = u_boundary(3,y)

# y=0 and y=4
for i in range(mx):
  x = 3.0*i/(mx-1)
  bc[0*mx + i] = u_boundary(x,0)
  bc[(my-1)*mx + i] = u_boundary(x,4)

# define initial values
initial_values = {}

for j in range(my):
  for i in range(mx):
    x = 3.0*i/(mx-1)
    y = 4.0*j/(my-1)
    initial_values[j*mx + i] = u_initial(x,y)

print(mx,my)

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "diffusion",                # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "Solvers": {
    "linearSolver":
    {
      "solverType":        "gamg",
      "preconditionerType": "none",
      "relativeTolerance": 1e-15,
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
      
      # gamg specific options:
      "gamgType": "classical",        # one of agg, geo, or classical 
      "cycleType": "cycleV",
      "nLevels": 25,
      
      "maxIterations": 1000,
      "dumpFormat": "default",
      "dumpFilename": "",
    }
  },

  "Heun" : {
    "initialValues": initial_values,
    "timeStepWidth": 1e-3,
    "endTime": 1.0,
    "timeStepOutputInterval": 100,
    "checkForNanInf": False,
    "inputMeshIsGlobal": True,
    "nAdditionalFieldVariables": 0,
    "additionalSlotNames": [],
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename": None,   # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    
    "FiniteElementMethod" : {
      # mesh parameters
      "nElements":        [nx, ny],
      "inputMeshIsGlobal": True,
      "physicalExtent":   [3.0, 4.0],
      "physicalOffset":   [0.0, 0.0],
      
      "solverName":       "linearSolver",
      "slotName": "",
      
      # problem parameters
      "prefactor":        0.5,
    },

    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 10, "filename": "out/filename", "binary": True, "fixedFormat": False, "onlyNodalValues": True, "combineFiles": True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "outputInterval": 10, "filename": "out/diffusion2d", "binary": True, "onlyNodalValues": True, "combineFiles": True, "fileNumbering": "incremental"}
    ]
  },
}
