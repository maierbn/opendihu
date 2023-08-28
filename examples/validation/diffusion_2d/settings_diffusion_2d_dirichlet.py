import numpy as np

# size of domain
size_x = 10
size_y = 10

# problem definition

# initial values
def u_initial(x, y):
  d = np.minimum(np.sqrt((x-2)**2 + (y-size_y/2)**2), np.pi)
  return 1 + np.cos(d)

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
  print(f"quadratic, mx={mx}")

# define Dirichlet boundary conditions
bc = {}
# x=0
for j in range(my):
  y = size_y*j/(my-1)
  bc[j*mx + 0] = 0

# define initial values
initial_values = {}

for j in range(my):
  for i in range(mx):
    x = size_x*i/(mx-1)
    y = size_y*j/(my-1)
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
    },
    "implicitSolver":
    {
      "solverType":        "lu",
      "preconditionerType": "none",
      "relativeTolerance": 1e-15,
      "absoluteTolerance": 1e-10
    }
  },

  "CrankNicolson" : {
    "initialValues": initial_values,
    "timeStepWidth": 1e-5,
    "endTime": 1.0,
    "timeStepOutputInterval": 100,
    "checkForNanInf": False,
    "inputMeshIsGlobal": True,
    "nAdditionalFieldVariables": 0,
    "additionalSlotNames": [],
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename": None,   # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    "solverName":       "implicitSolver",
    
    "FiniteElementMethod" : {
      # mesh parameters
      "nElements":        [nx, ny],
      "inputMeshIsGlobal": True,
      "physicalExtent":   [size_x, size_y],
      "physicalOffset":   [0.0, 0.0],
      
      "solverName":       "linearSolver",
      "slotName": "",
      
      # problem parameters
      "prefactor":        3,
    },

    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1000, "filename": "out_dirichlet/diffusion2d", "binary": True, "fixedFormat": False, "onlyNodalValues": True, "combineFiles": True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "outputInterval": 1000, "filename": "out_dirichlet/diffusion2d", "binary": True, "onlyNodalValues": True, "combineFiles": True, "fileNumbering": "incremental"}
    ]
  },
}
