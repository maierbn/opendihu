# Modified diffusion 1D 
# Compute 

n = 5   # number of elements  

config = {
  "logFormat": "csv",
  "Solvers": {
    "linearSolver": {
      "solverType": "gmres",          # the solver type, refer to PETSc documentation about implemented solvers and preconditioners (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html)
      "preconditionerType": "none",   # the preconditioner type (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType)
      "relativeTolerance": 1e-15,     # the relative tolerance of the linear solver
      "absoluteTolerance": 1e-10,     # 1e-10 absolute tolerance of the residual    
      "maxIterations": 1e4,           # maximum number of iterations of the linear solver
      "dumpFilename": "",             # a filename to dump the system matrix and rhs in every timestep, "" to disable
      "dumpFormat": "matlab",         # the file format of the dump file, "matlab", "ascii" or "default"
    }
  },
  "MyNewTimesteppingSolver": {        # this is the name of the solver, as given in the constructor to the timestepping object
    "myOption": 42,                   # example option that is parsed in the constructor 
    "option1": "blabla",              # another example option that is parsed in the data object
    
    # option for the timestepping of MyNewTimesteppingSolver
    "endTime": 10.0,                  # end time of the simulation
    "timeStepWidth": 2.0,             # time step width
    "timeStepOutputInterval": 1,      # how often to print the current timestep to console 
       
    # settings for the nested solver
    "ImplicitEuler" : {
      "numberTimeSteps": 5,
      "endTime": 0.1,
      "initialValues": [2,5,10,4,-2,2],    # the initial values
      "dirichletBoundaryConditions": {}, # Dirichlet boundary conditions as dict
      "dirichletOutputFilename":  None,  # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "inputMeshIsGlobal": True,         # initial values and BC's are given for all dofs, even if executed in parallel
      "timeStepOutputInterval": 1,       # how often to print the current timestep to console 
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
    },
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/paraview", "binary": True, "fixedFormat": False, "onlyNodalValues": True, "combineFiles": False},
      {"format": "PythonFile", "filename": "out/python", "outputInterval": 1, "binary":False, "onlyNodalValues": True}
    ]
  }
}
