MultidomainSolver
===================

This solves the multidomain equations given below. There is the `MultidomainSolver` which only solves the equations on the muscle domain 
(see example `examples/electrophysiology/multidomain/multidomain_no_fat`) 
and there is the `MultidomainWithFatSolver` which additionally considers a fat and skin layer. It uses a composite mesh of muscle and body domain.
(See example `examples/electrophysiology/multidomain/multidomain_with_fat`.)

The first and second Multidomain equation for compartments :math:`k = 1, \dots, N_\text{MU}` are given below:

  .. math::
    \color{red}{\textrm{div}\big(\sigma_e \,\textrm{grad}( \phi_e)\big) + \sum\limits_{k=1}^{N_\text{MU}} f_r^k\,\textrm{div}\big(\sigma_i^k\,\textrm{grad}(\phi_i^k)\big)  = 0}\\
    \color{red}{\textrm{div}\big(\sigma_i^k\,\textrm{grad}(\phi_i^k)\big)} = \color{orange}{ A_m^k\,\big(C_m^k \dfrac{\partial V_m^k}{\partial t} + I_\text{ion}(V_m^k, l_\text{HS}, \dot{l}_\text{HS}, \textbf{y}^k)\big),} \quad \forall k \in \{1, \dots, N_\text{MU}\}\\
    \color{orange}{\textbf{y}^k(t) = g(V_m^k, \textbf{y}^k(t))} \quad \forall k \in \{1, \dots, N_\text{MU}\}
  
  Reference: `Paper <https://link.springer.com/article/10.1007%2Fs10237-019-01214-5>`_
    
C++ code
----------

C++ code for a typical usage of the `MultidomainSolver`, i.e. without fat layer:

.. code-block:: c
  
  OperatorSplitting::Strang<
    Control::MultipleInstances<
      TimeSteppingScheme::Heun<
        CellmlAdapter<
          4,9,  // nStates,nAlgebraics: 57,1 = Shorten, 4,9 = Hodgkin Huxley
          FunctionSpace::FunctionSpace<MeshType,BasisFunction::LagrangeOfOrder<1>>  // same function space as for anisotropic diffusion
        >  
      >
    >,
    TimeSteppingScheme::MultidomainSolver<              // multidomain
      SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fibre directions
        MeshType,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<3>,
        Equation::Static::Laplace
      >,
      SpatialDiscretization::FiniteElementMethod<   // anisotropic diffusion
        MeshType,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<5>,
        Equation::Dynamic::DirectionalDiffusion
      >
    >
  >
  
C++ code for a typical usage of the `MultidomainWithFatSolver`, i.e. with fat layer:

.. code-block:: c
  
  OperatorSplitting::Strang<
    Control::MultipleInstances<
      TimeSteppingScheme::Heun<
        CellmlAdapter<
          4,9,  // nStates,nAlgebraics: 57,1 = Shorten, 4,9 = Hodgkin Huxley
          FunctionSpace::FunctionSpace<MeshType,BasisFunction::LagrangeOfOrder<1>>  // same function space as for anisotropic diffusion
        >  
      >
    >,
    TimeSteppingScheme::MultidomainWithFatSolver<       // multidomain
      SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fibre directions
        MeshType,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<3>,
        Equation::Static::Laplace
      >,
      SpatialDiscretization::FiniteElementMethod<       // anisotropic diffusion
        MeshType,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<5>,
        Equation::Dynamic::DirectionalDiffusion
      >,
      SpatialDiscretization::FiniteElementMethod<       // isotropic diffusion in fat layer
        MeshType,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<5>,
        Equation::Dynamic::IsotropicDiffusion
      >
    >
  >
  
Python settings
---------------------------

The MultidomainSolver is a time stepping scheme and therefor has the basic options of all time stepping schemes, such as ``timeStepWidth``, ``endTime``, ``timeStepOutputInterval``.

The options below can be used for both `MultidomainSolver` and `MultidomainWithFatSolver`. They are explained by the comments.

.. code-block:: python

  "MultidomainSolver": {
    "timeStepWidth":                    variables.dt_multidomain,             # time step width of the subcellular problem
    "endTime":                          variables.end_time,                   # end time, this is not relevant because it will be overridden by the splitting scheme
    "timeStepOutputInterval":           100,                                  # how often the output timestep should be printed
    "durationLogKey":                   "duration_multidomain",               # key for duration in log.csv file
    "slotNames":                        ["vm_old", "vm_new", "g_mu", "g_tot"],  # names of the data connector slots, maximum length per name is 10 characters. g_mu is gamma (active stress) of the compartment, g_tot is the total gamma
    
    # material parameters for the compartments
    "nCompartments":                    variables.n_compartments,             # number of compartments
    "compartmentRelativeFactors":       variables.relative_factors.tolist(),  # list of lists of (the factors for all dofs), because "inputIsGlobal": True, this contains the global dofs
    "inputIsGlobal":                    True,                                 # if values and dofs correspond to the global numbering
    "am":                               [variables.get_am(mu_no) for mu_no in range(variables.n_compartments)],   # Am parameter for every motor unit (ration of surface to volume of fibers)
    "cm":                               [variables.get_cm(mu_no) for mu_no in range(variables.n_compartments)],   # Cm parameter for every motor unit (capacitance of the cellular membrane)
    
    # solver options
    "solverName":                       "multidomainLinearSolver",            # reference to the solver used for the global linear system of the multidomain eq.
    "alternativeSolverName":            "multidomainAlternativeLinearSolver", # reference to the alternative solver, which is used when the normal solver diverges
    "subSolverType":                    "gamg",                               # sub solver when block jacobi preconditioner is used
    "subPreconditionerType":            "none",                               # sub preconditioner when block jacobi preconditioner is used
    #"subPreconditionerType":            "boomeramg",                          # sub preconditioner when block jacobi preconditioner is used, boomeramg is the AMG preconditioner of HYPRE

    # gamg specific options:
    "gamgType":                         "classical",                          # one of agg, geo, or classical 
    "cycleType":                        "cycleV",                             # either cycleV or cycleW
    "nLevels":                          25,
    "hypreOptions":                     "-pc_hypre_boomeramg_strong_threshold 0.7",       # additional options if a hypre preconditioner is selected
    
    "theta":                            variables.theta,                      # weighting factor of implicit term in Crank-Nicolson scheme, 0.5 gives the classic, 2nd-order Crank-Nicolson scheme, 1.0 gives implicit euler
    "useLumpedMassMatrix":              variables.use_lumped_mass_matrix,     # which formulation to use, the formulation with lumped mass matrix (True) is more stable but approximative, the other formulation (False) is exact but needs more iterations
    "useSymmetricPreconditionerMatrix": variables.use_symmetric_preconditioner_matrix,    # if the diagonal blocks of the system matrix should be used as preconditioner matrix
    "initialGuessNonzero":              variables.initial_guess_nonzero,      # if the initial guess for the 3D system should be set as the solution of the previous timestep, this only makes sense for iterative solvers
    "enableFatComputation":             True,                                 # disabling the computation of the fat layer is only for debugging and speeds up computation. If set to False, the respective matrix is set to the identity
    "showLinearSolverOutput":           variables.show_linear_solver_output,  # if convergence information of the linear solver in every timestep should be printed, this is a lot of output for fast computations
    "updateSystemMatrixEveryTimestep":  False,                                # if this multidomain solver will update the system matrix in every first timestep, use this only if the geometry changes, e.g. by contraction
    "updateSystemMatrixInterval":       1,                                    # if updateSystemMatrixEveryTimestep is True, how often the system matrix should be rebuild, in terms of calls to the solver. (E.g., 2 means every second time the solver is called)
    "recreateLinearSolverInterval":     0,                                    # how often the Petsc KSP object (linear solver) should be deleted and recreated. This is to remedy memory leaks in Petsc's implementation of some solvers. 0 means disabled.
    "rescaleRelativeFactors":           True,                                 # if all relative factors should be rescaled such that max Σf_r = 1
    "setDirichletBoundaryConditionPhiE":False,                                # (set to False) if the last dof of the extracellular space (variable phi_e) should have a 0 Dirichlet boundary condition. However, this makes the solver converge slower.
    "setDirichletBoundaryConditionPhiB":False,                                # (set to False) if the last dof of the fat layer (variable phi_b) should have a 0 Dirichlet boundary condition. However, this makes the solver converge slower.
    "resetToAverageZeroPhiE":           True,                                 # if a constant should be added to the phi_e part of the solution vector after every solve, such that the average is zero
    "resetToAverageZeroPhiB":           True,                                 # if a constant should be added to the phi_b part of the solution vector after every solve, such that the average is zero
    }
  
The list of values for `compartmentRelativeFactors` is for the symbol :math:`f_r^k` in the equations. The values should be the same on every rank. 
It is beneficial to compute the values once and store them in a cache file. Note that the number of nodes in total can be different if the same settings are used for different numbers of ranks. Therefore it is not easily possible to run the program serially, compute the cache of `compartmentRelativeFactors` and reuse it for all ranks.
Instead, the cache has to be created by a parallel run, but then only rank 0 should compute the values. This is done in the `helper.py` script of the examples `examples/electrophysiology/multidomain/multidomain_no_fat` and `examples/electrophysiology/multidomain/multidomain_with_fat`.

This means that using the MultidomainSolver is not so trivial. Therefore, the two given examples should be reused or copied if any new example is to be created. 
Not all settings will be explained again in the following as they are already given above.


`solverName` and `alternativeSolverName`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There is th problem with the `HYPRE` AMG solver that is sometimes converges very fast and sometimes diverges. To be able to use it nevertheless, the `alternativeSolverName` solver is automatically used after the main solver `solverName` diverged.

theta
^^^^^^^^^
Weighting factor of the implicit term in the Crank-Nicolson scheme:

.. math::

  \dfrac{u^{(t+1)} - u^{(t)}}{dt} = \theta \cdot rhs(u^{(t+1)},t+1) + (1-\theta) \cdot rhs(u^{(t)},t)

:math:`\theta=0.5` gives the classic, 2nd-order Crank-Nicolson scheme, :math:`\theta=1` leads to implicit euler. The fully implicit schemes was found to be more robust and should be used, i.e. ``theta=1``.

useLumpedMassMatrix
^^^^^^^^^^^^^^^^^^^^^^^^^
Which formulation to use, the formulation with lumped mass matrix (`True`) is more stable but approximative, the other formulation (`False`) is exact but needs more iterations. Usually, this option can be set to `True`.

useSymmetricPreconditionerMatrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the diagonal blocks of the system matrix should be used as preconditioner matrix. If set to false, the whole matrix is used for preconditioning.

initialGuessNonzero
^^^^^^^^^^^^^^^^^^^^^^^^^
If the initial guess for the 3D system is given by the solution of the previous timestep. This only makes sense for iterative solvers. A direct solver ``"lu"`` requires that this option is set to ``False``.

enableFatComputation
^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a switch to disable the fat layer in the multidomain formulation. Disabling the computation of the fat layer is only for debugging and speeds up computation. If set to False, the respective matrix is set to the identity.


showLinearSolverOutput
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If convergence information of the linear solver should be printedin every timestep. As this involves a lot of output for small and fast computations, it should be disabled. It can be useful for large and slow computations to see the, e.g., the number of iterations of the linear solver.


updateSystemMatrixEveryTimestep and updateSystemMatrixInterval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This option allows to create a new system matrix before every new solve. This is only required if the geometry changes, e.g., if a solid mechanics solver is deforming the domain.

If ``updateSystemMatrixEveryTimestep`` is set to `True`, the option ``updateSystemMatrixInterval`` determines, how frequently the system matrix will be rebuild. A value of 1 means in every first timestep of the multidomain solver, a higher value specifies every which call to the solver will compute a new system matrix. This is needed, e.g., if the multidomain solver and a muscle contraction solver are coupled and the multidomain solver is called more often than the muscle contraction solver.

recreateLinearSolverInterval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There appears to be a memory leak in some implementation of a PETSc solver that is visible during long runs. Using this option, it is possible to recreate the PETSc KSP object after the given number of time steps to free the memory. Apparently, the memory is still not freed despite deleting and recreating the PETSc solver.

rescaleRelativeFactors
^^^^^^^^^^^^^^^^^^^^^^^^^^^
If all relative factors should be rescaled such that max Σf_r = 1. 


setDirichletBoundaryConditionPhiE, setDirichletBoundaryConditionPhiB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the last dof of the extracellular space (variable phi_e) and the fat layer (variable phi_b) should have a 0 Dirichlet boundary condition, each.
This remove zero eigenvalues and makes the system regular.
However, this makes the solver converge slower. Therefore, these options are usually set to False.

resetToAverageZeroPhiE, resetToAverageZeroPhiB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If a constant should be added to the phi_e and phi_b parts of the solution vector after every solve, such that the average is zero.
This reduces the drift of the solution of the system is not regular. (If setDirichletBoundaryConditionPhi* is False).
