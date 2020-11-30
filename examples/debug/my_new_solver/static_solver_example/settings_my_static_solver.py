# Modified diffusion 1D 
# Compute 
import numpy as np

nx = 10
ny = 10

# boundary conditions (for linear elements)
# set bottom and top boundary to 1
bc = {}
for i in range(int(nx+1)):
  x = i/(nx+1.)
  #bc[i] = np.sin(x*np.pi)
  bc[i] = 1
  i2 = (nx+1)*ny + i
  #bc[i2] = np.sin(x*np.pi)
  bc[i2] = 1
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))

config = {
  "logFormat": "csv",
  "MyNewStaticSolver": {        # this is the name of the solver, as given in the constructor to the timestepping object
    "myOption": 42,                   # example option that is parsed in the constructor 
    "option1": "blabla",              # another example option that is parsed in the data object
    
    # settings for the nested solver
    "FiniteElementMethod" : {
      # mesh
      "nElements": [nx, ny],               # number of elements in x and y direction
      "physicalExtent": [1.0, 1.0],        # the physical size of the domain
      
      # solver
      "solverType": "gmres",          # the solver type, refer to PETSc documentation about implemented solvers and preconditioners (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html)
      "preconditionerType": "none",   # the preconditioner type (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType)
      "relativeTolerance": 1e-15,     # the relative tolerance of the linear solver
      "absoluteTolerance": 1e-10,     # 1e-10 absolute tolerance of the residual    
      "maxIterations": 1e4,           # maximum number of iterations of the linear solver
      "dumpFilename": "out/",             # a filename to dump the system matrix and rhs in every timestep, "" to disable
      "dumpFormat": "matlab",         # the file format of the dump file, "matlab", "ascii" or "default"
      
      # problem definition
      "prefactor": 1,                      # prefactor of the laplace operator
      "dirichletBoundaryConditions": bc,   # Dirichlet boundary conditions as dict  
      "dirichletOutputFilename":  None,    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "neumannBoundaryConditions": [{"element": j*nx, "constantVector": 1.0, "face": "0-"} for j in range(ny)],   # Neumann boundary conditions: at boundary
      "inputMeshIsGlobal": True,           # boundary condition values are given for all dofs, even if executed in parallel
    },

    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/paraview", "binary": True, "fixedFormat": False, "onlyNodalValues": True, "combineFiles": False},
      {"format": "PythonFile", "filename": "out/python", "outputInterval": 1, "binary":False, "onlyNodalValues": True}
    ]
  }
}
