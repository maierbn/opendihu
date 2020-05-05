# Diffusion 1D

n = 17   # number of elements

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "Solvers": {
    "linearSolver": {
      "solverType": "gmres",          # the solver type, refer to PETSc documentation about implemented solvers and preconditioners (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html)
      "preconditionerType": "none",   # the preconditioner type (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType)
      "relativeTolerance": 1e-15,     # the relative tolerance of the linear solver
      "absoluteTolerance": 1e-10,     # 1e-10 absolute tolerance of the residual norm
      "maxIterations": 1e4,           # maximum number of iterations of the linear solver
      "dumpFilename": "",             # a filename to dump the system matrix and rhs in every timestep, "" to disable
      "dumpFormat": "matlab",         # the file format of the dump file, "matlab", "ascii" or "default"
    }
  },
  "ImplicitEuler" : {
     "numberTimeSteps": 10,
     "startTime": 0.0,
     "endTime": 100,
     "initialValues": [2,2,4,5,2,2,2,0,0],    # the initial values
     "dirichletBoundaryConditions": {}, # Dirichlet boundary conditions as dict
     "inputMeshIsGlobal": True,         # initial values and BC's are given for all dofs, even if executed in parallel
     "timeStepOutputInterval": 100000,       # how often to print the current timestep to console
     "nAdditionalFieldVariables": 0,    # for more complex nested solvers, the number of additional field variables that will be transferred without being touched
     "solverName": "linearSolver",      # the solver to use, referes to what was defined under "Solvers"

     "FiniteElementMethod" : {
        # mesh
        "nElements": n,                 # number of elements
        "physicalExtent": 4.0,          # the physical size of the domain

        "solverName": "linearSolver",   # the solver to use, referes to what was defined under "Solvers"
        "prefactor": 5.0,               # the prefactor 'c' of 'du/dt = c du^2/dx^2'
        "inputMeshIsGlobal": True,      # boundary conditions are given as global indices
        "nAdditionalFieldVariables": 0, # for more complex nested solvers, the number of additional field variables that will be transferred without being touched
     },
     "OutputWriter" : [
       #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "onlyNodalValues": True},
       {"format": "PythonFile", "filename": "out/diffusion1d_implicit", "outputInterval": 1, "binary":False, "onlyNodalValues": True}
     ]
  }
}
