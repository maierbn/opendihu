import numpy as np
import sys

# Problem definition
# ------------------
# 3D Poisson equation -Δu + f = 0
# Domain Ω = [0,2]x[0,3]x[0,4]
# Dirichlet boundary conditions
# u(a) = ua, u(b) = ub

# parameters
n = 1  # n elements

# if an extra command line argument is given, set number of elements to this value
if len(sys.argv) == 3:
  n = int(sys.argv[0])

# boundary conditions
def u(x,y,z):
  return x**3*y**2*z + 4*x**2*y**2*z**3 + 2*x*y**3*z - y*z**2 + 3*x**2*y + 1

def f(x,y,z):
  return 2*x**3*z + 24*x**2*y**2*z + 8*x**2*z**3 + 6*x*y**2*z + 12*x*y*z + 8*y**2*z**3 + 4*y


# discretization of problem
# -------------------------
import opendihu
if "regular" in opendihu.program_name:
  nx = 2*n
  ny = 3*n
  nz = 4*n

elif "structured" in opendihu.program_name:
  nx = n
  ny = n
  nz = n

else:
  print(f"Invalid program name {opendihu.program_name}")

if "linear" in opendihu.program_name:

  # number of nodes
  mx = nx + 1
  my = ny + 1
  mz = nz + 1

elif "quadratic" in opendihu.program_name:

  # number of nodes
  mx = 2*nx + 1
  my = 2*ny + 1
  mz = 2*nz + 1

else:
  print(f"Invalid program name {opendihu.program_name}")

# boundary conditions
bc = {}

# x=0 and x=2
for k in range(mz):
  for j in range(my):
    y = 3*j/(my-1)
    z = 4*k/(mz-1)
    bc[k*mx*my + j*mx + 0] = u(0,y,z)
    bc[k*mx*my + j*mx + mx-1] = u(2,y,z)

# y=0 and y=3
for k in range(mz):
  for i in range(mx):
    x = 2*i/(mx-1)
    z = 4*k/(mz-1)
    bc[k*mx*my + 0*mx + i] = u(x,0,z)
    bc[k*mx*my + (my-1)*mx + i] = u(x,3,z)

# z=0 and z=4
for j in range(my):
  for i in range(mx):
    x = 2*i/(mx-1)
    y = 3*j/(my-1)
    bc[0*mx*my + j*mx + i] = u(x,y,0)
    bc[(mz-1)*mx*my + j*mx + i] = u(x,y,4)

# right hand side
rhs = np.zeros(mz*my*mx)
for k in range(mz):
  for j in range(my):
    for i in range(mx):
      x = 2*i/(mx-1)
      y = 3*j/(my-1)
      z = 4*k/(mz-1)
      
      rhs[k*mx*my + j*mx + i] = f(x,y,z)

print(f"{mx} {my} {mz}")

output_filename = opendihu.program_name.replace(".py","")

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "poisson",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
    # mesh parameters
    "nElements": [nx, ny, nz],
    "physicalExtent": [2, 3, 4],
    "physicalOffset": [0, 0, 0],
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
      {"format": "Paraview", "interval": 1, "filename": "out/p", "binary": True, "fixedFormat": False},
      {"format": "PythonFile", "filename": f"out/{output_filename}", "outputInterval": 1, "binary": True, "onlyNodalValues": False, "fileNumbering": "incremental"}
    ]
  },
}
