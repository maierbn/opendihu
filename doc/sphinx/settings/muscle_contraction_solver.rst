Muscle contraction solver
===========================

This class solves a dynamic or quasi-static hyperelasticity problem, possibly with active stress terms. 
It adds connector slots, e.g. for the fiber stretch, :math:`\lambda`. It also contains an `OutputWriter` that writes the most important field variables to a single output file.
Internally it encapsulates the classes `HyperelasticitySolver` and `DynamicHyperelasticitySolver` with `Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D` as Material.
Currently, the Material is hardcoded, because it is specific to the active stress term. It should not be too difficult to make it also a template argument in the C++ core code.

The differences to `DynamicHyperelasticitySolver` is the following:

* Allows to choose dynamic or quasi-static problem (by option ``dynamic``).
* Adds the active stress values to the 2nd Piola-Kirchhoff stresses. This value is provided over the connector slot "γ". The active stress is computed as 

  .. math::
    
    \dfrac{1}{\lambda}  P_\text{max} f(\lambda/\lambda_\text{opt}) \gamma
    
  with the force-length relation

  .. math::

    f(\lambda/\lambda_\text{opt}) = -\dfrac{25}{4} \left(\dfrac{\lambda}{\lambda_\text{opt}}\right)^2 + \dfrac{25}{2} \dfrac{\lambda}{\lambda_\text{opt}} - 5.25,\\
    \lambda_\text{opt} = 1.2
    
* Computes fiber stretch, :math:`\lambda`, contraction velocity, :math:`\dot\lambda` (not yet implemented), and material traction, :math:`T`, and provides these values via a connector slot.
  The order of the slots is: λ, λdot, γ, T
* Provides an output writer which contains all important information. The `HyperelasticitySolver` and `DynamicHyperelasticitySolver` also provide OutputWriters,  but several with different data. Here, all static and dynamic quantities of interest, such as stresses, displacements and velocities are combined into a single file.
  
C++ instantiation
-----------------

.. code-block:: c

  // either with composite mesh:
  MuscleContractionSolver<
    Mesh::CompositeOfDimension<3>
  >
  
  // or with structured mesh:
  MuscleContractionSolver<
    Mesh::StructuredDeformableOfDimension<3>
  >

Alternatively, the `Term`, i.e. material can be specified as second template argument:

.. code-block:: c

  // either with composite mesh:
  MuscleContractionSolver<
    Mesh::CompositeOfDimension<3>,
    Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D
  >
    
The given value of `Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D` is the default but can be replaced. 
The material has to set ``usesActiveStress=true`` in order for the active stress to work.

And additionally, the parameter `withLargeOutputFiles` can be specified as the third parameter (default is true). If set to false, smaller output files that don't contain :math:`P,\dot{F}` and :math:`F` values will be created.

.. code-block:: c

  // either with composite mesh:
  MuscleContractionSolver<
    Mesh::CompositeOfDimension<3>,
    Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D,
    true
  >

Python settings
-----------------

In the following all possible options are listed. and explained by the comments.
Depending on the value of ``dynamic``, either an instance of ``DynamicHyperelasticitSolver`` or ``HyperlasticitySolver`` is included. See :doc:`hyperelasticity` and :doc:`dynamic_hyperelasticity` for details.

.. code-block:: python
  
  "MuscleContractionSolver": {
    "numberTimeSteps":              1,                         # number of timesteps to use per call
    "timeStepOutputInterval":       100,                       # How often the current timestep will be displayed, if this is >100 and numberTimeSteps is 1, nothing will be printed
    "Pmax":                         variables.pmax,            # maximum PK2 active stress
    "enableForceLengthRelation":    True,                      # if the factor f_l(λ_f) modeling the force-length relation (as in Heidlauf2013) should be multiplied. Set to false if this relation is already considered in the CellML model.
    "lambdaDotScalingFactor":       1.0,                       # scaling factor for the output of the lambda dot slot, i.e. the contraction velocity. Use this to scale the unit-less quantity to, e.g., micrometers per millisecond for the subcellular model.
    "OutputWriter" : [                                         # This is an output writer that writes files with all required fields.
      {"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D), "filename": "out/" + variables.scenario_name + "/mechanics_3D", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    ],
    "mapGeometryToMeshes":          [],                        # the mesh names of the meshes that will get the geometry transferred
    "slotNames":                    ["lambda", "ldot", "gamma", "T"],    # names of the connector slots, maximum 10 characters per name 
    "dynamic":                      True,                      # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem
    
    # the actual solid mechanics solver, this is either "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on the value of "dynamic"
    "DynamicHyperelasticitySolver": {
      "timeStepWidth":              variables.dt_3D,           # time step width 
      "durationLogKey":             "nonlinear",               # key to find duration of this solver in the log file
      "timeStepOutputInterval":     1,                         # how often the current time step should be printed to console
      
      "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
      "density":                    variables.rho,             # density of the material
      "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
      "residualNormLogFilename":    "log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
      "useAnalyticJacobian":        variables.use_analytic_jacobian,      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
      "useNumericJacobian":         not variables.use_analytic_jacobian,  # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
        
      "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
      # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
      
      # mesh
      "inputMeshIsGlobal":          True,                     # the mesh is given locally
      "meshName":                   "3Dmesh_quadratic",        # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
      "fiberMeshNames":             variables.fiber_mesh_names,  # fiber meshes that will be used to determine the fiber direction, for multidomain there are no fibers so this would be empty list
      #"fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system

      # solving
      "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
      #"loadFactors":                [0.25, 0.66, 1.0],                # load factors for every timestep
      "loadFactors":                [],                        # no load factors, solve problem directly
      "loadFactorGiveUpThreshold":   1,                      # when to abort the solve
      "scaleInitialGuess":          False,                      # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
      "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
      
      # boundary and initial conditions
      "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
      "neumannBoundaryConditions":   variables.elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
      "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
      "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
      "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
      
      "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
      "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
      "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
      "constantBodyForce":           variables.constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
      
      "dirichletOutputFilename":    "out/"+variables.scenario_name+"/dirichlet_boundary_conditions",     # output filename for the dirichlet boundary conditions, set to "" to have no output
      "totalForceLogFilename":       "out/"+variables.scenario_name+"/tendon_force.csv",              # filename of a log file that will contain the total (bearing) forces and moments at the top and bottom of the volume
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
          #{"format": "Paraview", "outputInterval": int(output_interval/dt), "filename": "out/dynamic", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          #{"format": "Paraview", "outputInterval": int(1./variables.dt_3D*variables.output_timestep_3D), "filename": "out/"+variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        ],
      },
      # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
      "LoadIncrements": {   
        "OutputWriter" : [
          #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
        ]
      },
    }
  }
  
See the description in :doc:`dynamic_hyperelasticity` for all details, such as boundary conditions.


timeStepOutputInterval
^^^^^^^^^^^^^^^^^^^^^^^^^^^
How often the current timestep will be printed to the console and log. 
If this is >100 and numberTimeSteps is 1, nothing will be printed.


Pmax
^^^^^^^
The maximum active stress value, also called :math:`S_\text{max}` in the derivation.
    
enableForceLengthRelation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the function :math:`f_l(λ_f)` which models the force-length relation (as in Heidlauf2013) should be multiplied to the active stress value.
Set to false if this relation is already considered in the CellML model.

lambdaDotScalingFactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A scaling factor for the output of the lambda dot slot, i.e. the contraction velocity. 
Use this to scale the unit-less quantity to, e.g., micrometers per millisecond for the subcellular model. Usually this factor can be set to 1.

mapGeometryToMeshes
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a list of opendihu mesh names that will be updated with the deformed geometry of the 3D mesh.
If this list is not empty, only the specified meshes will be deformed according to the computed displacements.

If this list is empty, the meshes of all connected slots will automatically be deformed, as this is usually what you want.

slotNames
^^^^^^^^^^^^^^
A list names of the connector slots, maximum 6 characters per name, see :doc:`output_connector_slots` for details.
    
dynamic
^^^^^^^^^^^
This is a bool variable that specifies if the static solver (:doc:`hyperelasticity`) or the dynamic solver (:doc:`dynamic_hyperelasticity`) should be used.
   
