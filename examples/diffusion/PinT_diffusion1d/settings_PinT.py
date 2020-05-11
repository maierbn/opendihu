# Modified diffusion 1D
# Compute

import numpy as np

n = 9   # number of elements

config = {
"solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
 "Solvers": {
   "linearSolver".format(k): {
     "solverType": "gmres",          # the solver type, refer to PETSc documentation about implemented solvers and preconditioners (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html)
     "preconditionerType": "none",   # the preconditioner type (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType)
     "relativeTolerance": 1e-15,     # the relative tolerance of the linear solver
     "absoluteTolerance": 1e-10,     # absolute tolerance of the linear solver
     "maxIterations": 1e4,           # maximum number of iterations of the linear solver
     "dumpFilename": "",             # a filename to dump the system matrix and rhs in every timestep, "" to disable
     "dumpFormat": "matlab",         # the file format of the dump file, "matlab", "ascii" or "default"
   } for k in range(n)
  },
  # for k in range (n):
  "Meshes":{
          #"mesh_0": {
          #  "nElements": 2,                 # number of elements
          #  "physicalExtent": 4.0,          # the physical size of the domain
          #  "inputMeshIsGlobal": True,      # boundary conditions are given as global indices
          #},
           "mesh_{}".format(k): {
              "nElements": [(2 ** (k))] if k > 0 else [2],                 # number of elements
              # "nElements": 32,                 # number of elements
              "physicalExtent": 4.0,          # the physical size of the domain
              "inputMeshIsGlobal": True,      # boundary conditions are given as global indices
            } for k in range(n)
        },
  "PinTIE": {        # this is the name of the solver, as given in the constructor to the timestepping object
    "tstart": 0,                    # Start time
    "tstop": 100,                     # End time
    "ntime": 1000000,                      # number of time steps
    "nspace":   32,
    "Initial Guess": [2,2,4,5,2,2,2,0, 2,2,4,5,2,2,2,0,2,2,4,5,2,2,2,0,2,2,4,5,2,2,2,0],
    "option1": "blabla",              # another example option that is parsed in the data object
    "nRanksInSpace": 1,            # number of processes that compute the spatial domain in parallel
    "TimeSteppingScheme": [
    {
      "ImplicitEuler": {
        "numberTimeSteps": 1,
        "startTime":0.0,
        "endTime": j,
        # "timeStepWidth": 1,
        "initialValues": [2,2,4,5,2,2,2,0, 2,2,4,5,2,2,2,0,2,2,4,5,2,2,2,0,2,2,4,5,2,2,2,0],    # the initial values
        "dirichletBoundaryConditions": {}, # Dirichlet boundary conditions as dict
        "inputMeshIsGlobal": True,         # initial values and BC's are given for all dofs, even if executed in parallel
        "timeStepOutputInterval": 1,       # how often to print the current timestep to console
        "nAdditionalFieldVariables": 0,    # for more complex nested solvers, the number of additional field variables that will be transferred without being touched
        "solverName": "linearSolver".format(j),      # the solver to use, referes to what was defined under "Solvers"
        "FiniteElementMethod" : {
           "meshName": "mesh_{}".format(j),
           "solverName": "linearSolver".format(j),   # the solver to use, referes to what was defined under "Solvers"
           "prefactor": 5.0,               # the prefactor 'c' of 'du/dt = c du^2/dx^2'
           "inputMeshIsGlobal": True,      # boundary conditions are given as global indices
           "nAdditionalFieldVariables": 0, # for more complex nested solvers, the number of additional field variables that will be transferred without being touched
        },
        "OutputWriter" : [
          # {"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False, "onlyNodalValues": True},
          {"format": "PythonFile", "filename": "out/diffusion1d_PinT", "outputInterval": 1, "binary":False, "onlyNodalValues": True}
        ]if j == 20 else [],
      }
    } for j in range(n)],
  },
 }
