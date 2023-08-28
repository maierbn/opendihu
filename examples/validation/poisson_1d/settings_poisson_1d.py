import numpy as np
import sys

# Problem definition
# ------------------
# 1D Poisson equation -Δu + f = 0
# Domain Ω = [a, b] with a=0
# Dirichlet boundary conditions
# u(a) = ua, u(b) = ub

# parameters
n = 11  # n elements

# if an extra command line argument is given, set number of elements to this value
if len(sys.argv) == 3:
  n = int(sys.argv[0])

# boundary conditions
b = 3
ua = 1 
ub = 2

# right hand side
def f(x):
  return 1-x**2

def F(x):
  return x-1/3*x**3 + 13/(4*b)

def fprime(x):
  return -2*x

# discretization of problem
# -------------------------
import opendihu
if "linear" in opendihu.program_name:

  # number of nodes
  m = n+1 

  # boundary conditions
  bc = {0: ua, m-1: ub}

  # right hand side
  rhs = [f(b*i/(m-1)) for i in range(m)]

elif "quadratic" in opendihu.program_name:

  # number of nodes
  m = 2*n+1

  # boundary conditions
  bc = {0: ua, m-1: ub}

  # right hand side
  rhs = [f(b*i/(m-1)) for i in range(m)]

elif "hermite" in opendihu.program_name:
  
  # for Hermite discretization, every node has to dofs: the function value and the derivative
  m = 2*n+2     # number of nodes
  mx = b/n      # element size

  # boundary conditions
  bc = {0: ua, 1: F(0)*mx, m-2: ub}  # Dirichlet bc
  
  # set right hand side (coefficients to the ansatz functions)
  rhs = [f(b*i/(m-2)) if i%2 == 0 else fprime(b*(i-2)/(m-2)) for i in range(m)]

output_filename = opendihu.program_name.replace(".py","")

print(f"m={m}, bc: {bc}, rhs: {rhs}")

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "poisson",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
    # mesh parameters
    "nElements": [n],
    "physicalExtent": [b],
    "inputMeshIsGlobal": True,
    
    # problem parameters
    "prefactor": 1,
    "rightHandSide": list(rhs),         # provide the rhs vector as list, not as numpy array
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename": None,    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    
    # solver parameters
    "solverType": "lu",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
    "maxIterations": 1000,
    "dumpFormat": "ascii",
    "dumpFilename": "out/poisson",
    "slotName": "",
    
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out/p", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": f"out/{output_filename}", "outputInterval": 1, "binary": True, "onlyNodalValues": False, "fileNumbering": "incremental"}
    ]
  },
}
