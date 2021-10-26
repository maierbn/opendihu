# Diffusion 1D

n = 5   # number of elements  

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "diffusion",                # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
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
     "numberTimeSteps": 5,
     "endTime": 0.1,
     "initialValues": [2,2,4,5,2,2],    # the initial values
     "dirichletBoundaryConditions": {}, # Dirichlet boundary conditions as dict
     "dirichletOutputFilename": None,   # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
     "inputMeshIsGlobal": True,         # initial values and BC's are given for all dofs, even if executed in parallel
     "timeStepOutputInterval": 1,       # how often to print the current timestep to console 
     "nAdditionalFieldVariables": 0,    # for more complex nested solvers, the number of additional field variables that will be transferred without being touched
     "additionalSlotNames": [],         # the slot names for the additional field variables, this is a list of strings with maximum 10 characters each
     "solverName": "linearSolver",      # the solver to use, referes to what was defined under "Solvers"
     "checkForNanInf": True,            # if the solution should be checked for nan or inf values which indicate an instability, this may be an expensive check
        
     "FiniteElementMethod" : {
        # mesh
        "nElements": n,                 # number of elements
        "physicalExtent": 4.0,          # the physical size of the domain
        
        "solverName": "linearSolver",   # the solver to use, referes to what was defined under "Solvers"
        "prefactor": 5.0,               # the prefactor 'c' of 'du/dt = c du^2/dx^2'
        "inputMeshIsGlobal": True,      # boundary conditions are given as global indices
        "slotName": "",
     },
     "OutputWriter" : [
       #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "onlyNodalValues": True, "fileNumbering": "incremental"},
       {"format": "PythonFile", "filename": "out/diffusion1d_implicit", "outputInterval": 1, "binary":False, "onlyNodalValues": True, "fileNumbering": "incremental"}
     ]
  }
}
