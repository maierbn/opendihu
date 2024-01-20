Dynamic hyperelasticity
===========================

The dynamic hyperelasticity solver solves nonlinear, hyperelastic, dynamic problems of solid mechanics in 3D. Internally, the class inherits from the :doc:`hyperelasticity` solver.
All python settings of the hyperelasticity solver as well as the method to specify the material is shared between the dynamic and the static hyperelasticity solvers.

C++ instantiation
-----------------

.. code-block:: c

  SpatialDiscretization::DynamicHyperelasticitySolver<
    Material
  >
  // or:
  SpatialDiscretization::DynamicHyperelasticitySolver<Material, true>   // default, same as without "false"
  SpatialDiscretization::DynamicHyperelasticitySolver<Material, false>

The two template parameters are same as in :doc:`hyperelasticity`. The first, ``Material`` is a class that describes the used constitutive equations at compile time.
The second specifies if there should be also the :math:`P` and :math:`F` field variables in the output files, which produces larger files.

Python settings
-----------------

In the following all possible options for the dynamic hyperelasticity solver are listed. They are explained by the comments. 
Most of them are also present in the :doc:`HyperelasticitySolver <hyperelasticity>`. Information regarding Neumann and Dirichlet boundary conditions, which correspond to traction/forces and prescribed nodes can be found under :doc:`boundary_conditions`.

.. code-block:: python

  "DynamicHyperelasticitySolver": {
    "timeStepWidth":              variables.dt_elasticity,      # time step width 
    "endTime":                    variables.end_time,           # end time of the simulation time span    
    "durationLogKey":             "duration_mechanics",         # key to find duration of this solver in the log file
    "timeStepOutputInterval":     1,                            # how often the current time step should be printed to console
    
    "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
    "density":                    variables.rho,                # density of the material
    "displacementsScalingFactor": 1.0,                          # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
    "residualNormLogFilename":    "log_residual_norm.txt",      # log file where residual norm values of the nonlinear solver will be written
    "useAnalyticJacobian":        True,                         # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
    "useNumericJacobian":         False,                        # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables":   False,                        # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "meshName":                   "3Dmesh_quadratic",           # mesh with quadratic Lagrange ansatz functions
    "inputMeshIsGlobal":          True,                         # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
  
    "fiberMeshNames":             [],                           # fiber meshes that will be used to determine the fiber direction
    "fiberDirection":             [],                           # if fiberMeshNames is empty, directly set the constant fiber direction, in global coordinate system
    "fiberDirectionInElement":    [0,0,1],                      # if fiberMeshNames and fiberDirections are empty, directly set the constant fiber direction, in element coordinate system
    
    # nonlinear solver
    "relativeTolerance":          1e-5,                         # 1e-10 relative tolerance of the linear solver
    "absoluteTolerance":          1e-10,                        # 1e-10 absolute tolerance of the residual of the linear solver       
    "solverType":                 "preonly",                    # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType":         "lu",                         # type of the preconditioner
    "maxIterations":              1e4,                          # maximum number of iterations in the linear solver
    "snesMaxFunctionEvaluations": 1e8,                          # maximum number of function iterations
    "snesMaxIterations":          10,                           # maximum number of iterations in the nonlinear solver
    "snesRelativeTolerance":      1e-5,                         # relative tolerance of the nonlinear solver
    "snesLineSearchType":         "l2",                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesAbsoluteTolerance":      1e-5,                         # absolute tolerance of the nonlinear solver
    "snesRebuildJacobianFrequency": 5,                          # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
    
    #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),          # dump system matrix and right hand side after every solve
    "dumpFilename":               "",                           # dump disabled
    "dumpFormat":                 "matlab",                     # default, ascii, matlab
    
    #"loadFactors":                [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
    #"loadFactors":                [0.5, 1.0],                   # load factors for every timestep
    "loadFactors":                [],                           # no load factors, solve problem directly
    "loadFactorGiveUpThreshold":  4e-2,                         # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the progression between two subsequent load factors gets smaller than this value, the solution is aborted.
    "scaleInitialGuess":          False,                        # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
    "nNonlinearSolveCalls":       1,                            # how often the nonlinear solve should be called
    
    # boundary and initial conditions
    "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
    "neumannBoundaryConditions":   variables.elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
    "updateNeumannBoundaryConditionsFunction": None,                    # function that updates the Neumann BCs while the simulation is running
    "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step
    
    "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
    "constantBodyForce":           variables.constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
    
    "dirichletOutputFilename":     "out/"+scenario_name+"/dirichlet_boundary_conditions_tendon",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "totalForceLogFilename":       "out/muscle_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
    "totalForceLogOutputInterval": 10,                                  # output interval when to write the totalForceLog file
    "totalForceBottomElementNosGlobal":  [j*nx + i for j in range(ny) for i in range(nx)],                  # global element nos of the bottom elements used to compute the total forces in the log file totalForceLogFilename
    "totalForceTopElementNosGlobal":     [(nz-1)*ny*nx + j*nx + i for j in range(ny) for i in range(nx)],   # global element nos of the top elements used to compute the total forces in the log file totalForceTopElementsGlobal

    
    # define which file formats should be written
    # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
    "OutputWriter" : [
      
      # Paraview files
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      
      # Python callback function "postprocess"
      #{"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": ""},
    ],
    # 2. additional output writer that writes also the hydrostatic pressure
    "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
      "OutputWriter" : [
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
    # 3. additional output writer that writes virtual work terms
    "dynamic": {    # output of the dynamic solver, has additional virtual work values 
      "OutputWriter" : [   # output files for displacements function space (quadratic elements)
        {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/dynamic", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ],
    },
    # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  }
  
  
The following options only apply to ``DynamicHyperelasticitySolver`` and not ``HyperelasticitySolver``:

`timeStepWidth`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The time step width of the time stepping solver, i.e., the timesteps in which the elasticity problem gets solved.

`endTime`
^^^^^^^^^^^^^^^^^^^
End time of the simulation.

`timeStepOutputInterval`
^^^^^^^^^^^^^^^^^^^^^^^^^
In which interval the current timestep will be written to the console.

`density`
^^^^^^^^^^^^
A constant density of the body, needed for the inertia effects.


  
`updateDirichletBoundaryConditionsFunction` (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a callback function that will be called regularly, in the interval given by the parameter `updateDirichletBoundaryConditionsFunctionCallInterval`. 
It allows to set new values for the Dirichlet boundary conditions, i.e. prescribed displacements and velocities. Set it to `None` to disable the callback.

The call function has the following form:

.. code-block:: python

  dirichlet_bc = {}
  dirichlet_bc[2] = 1   # prescribed dof 2 to value 1

  # Function to update dirichlet boundary conditions over time, t.
  # This function returns "dirichlet_bc". Only those entries can be updated that were also initially set.
  def update_dirichlet_boundary_conditions(t):
    
    dirichlet_bc[2] = 4   # change prescribed value of dof 2 to be value 4
    return dirichlet_bc

The only given argument, ``t``, is the current simulation time. The return value has to be a dict in the format that fits the parameter `dirichletBoundaryConditions`.
It is recommended to use a global variable, e.g. named ``dirichlet_bc``, that holds such a dict with all Dirichlet boundary conditions. 
Then, in the callback function, this variable is modified and returned.

Only the entries which were initially set can be modified. The reason for this is, that the prescribed dofs affect the matrix structure and the system matrix will not be reformed every time this callback was called, because this would be too expensive.

`updateDirichletBoundaryConditionsFunctionCallInterval` (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This option is the interval in which the callback function `updateDirichletBoundaryConditionsFunction` will be called. Only if `updateDirichletBoundaryConditionsFunction` was given in the config, this option is mandatory.


`updateNeumannBoundaryConditionsFunction` (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a callback function that will be called regularly, in the interval given by the parameter `updateDirichletBoundaryConditionsFunctionCallInterval`. 
It allows to set new Neumann boundary conditions, i.e. surface traction values. Set it to `None` to disable the callback.

The callback function has the following form:

.. code-block:: python

  # Function to update Neumann boundary conditions over time
  def update_neumann_boundary_conditions(t):
    
    # Neumann boundary conditions
    k = 0
    factor = np.sin(t/10. * 2*np.pi)*0.1
    neumann_bc = [{"element": k*nx*ny + j*nx + i, "constantVector": [factor,0,0], "face": "2-", "isInReferenceConfiguration": True} for j in range(ny) for i in range(nx)]
    #neumann_bc = []

    config = {
      "inputMeshIsGlobal": True,
      "divideNeumannBoundaryConditionValuesByTotalArea": False,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
      "neumannBoundaryConditions": neumann_bc
    }
    
    print("update neumann bc for t={}: {}".format(t,config))
    return config

The only given argument, ``t``, is the current simulation time. The return value has to be a config dict in the format shown above.
The options `"inputMeshIsGlobal"`, `"divideNeumannBoundaryConditionValuesByTotalArea"` and `"neumannBoundaryConditions"` have the same meaning as in the normal `config`.
 
This means the value of "neumannBoundaryConditions" has the usual list format for Neumann boundary conditions.

Changing Neumann boundary condition values only affects the right hand side of the mechanics problem. Therefore, Neumann BC values can be set for any number of elements, unlike in the Dirichlet BC callback. Previous Neumann boundary conditions are deleted.

`updateNeumannBoundaryConditionsFunctionCallInterval` (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This option is the interval in which the callback function `updateNeumannBoundaryConditionsFunction` will be called.
Only if `updateNeumannBoundaryConditionsFunction` was given in the config, this option is mandatory.


