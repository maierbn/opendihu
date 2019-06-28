
Solver
=======

A linear solver is e.g. needed for FiniteElementMethod class, for MultidomainSolver etc. whenever a linear system needs to be solved. It usually consists of a preconditioner and the actual solver. The preconditioners and solvers are provided by PETSc.

In the python settings, solvers can be defined at the beginning of the ``config`` dictionary under the ``"Solvers"`` key, similar to meshes under the ``"Meshes"`` key.
The item ``"Solvers"`` is itself a dictionary, where the properties of every solver are listed with a *solver name* as key. The solver name can be chosen abitrarily.
It is used to reference the solver later, where it is needed, e.g. in a FiniteElementMethod object.

.. code-block:: python

  config = {
    "Solvers": {
      "activationSolver": {
        "solverType": "gmres",
        "preconditionerType": "none",
        "relativeTolerance": 1e-5,
        "maxIterations": 1e4,
      },
      "otherSolver": {
         # properties of this solver
      }
    }
  }
  
solverType
~~~~~~~~~~~
*Default: gmres*

The *KSPType* of the solver, i.e. which solver to use. Possible values are:

- gmres
- sor
- jacobi
- cg
- bcgs
- lu
- cholesky
- gamg
- richardson
- chebyshev
- preonly (This means only use the preconditioner)

See `this PETSc page <https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType>`_ for details. All strings defined there are also possible.

preconditionerType
~~~~~~~~~~~~~~~~~~~
*Default: none*

The preconditioner type of PETSc to use. Possible values are:

- jacobi
- sor
- lu
- ilu  (incomplete LU factorization)
- gamg (geometric algebraic multigrid)
- none

See `the PETSc page on PCType <https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html>`_ for more information. All strings defined there are also possible.

relativeTolerance
~~~~~~~~~~~~~~~~~~
The relative tolerance of the residuum after which the solver is converged. Relative means relative to the initial residuum. 

maxIterations
~~~~~~~~~~~~~~
The maximum number of iterations after which the solver aborts and states divergence.

Internally the solver is configured with

.. code-block:: c
  
  //                        relative tol,      absolute tol,  diverg tol.,   max_iterations
  KSPSetTolerances (*ksp_, relativeTolerance_, PETSC_DEFAULT, PETSC_DEFAULT, maxIterations_)

See the `PETSc documentation for KSPSetTolerances <https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html>`_ to understand what that means.

Command line options for PETSc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PETSc uses its own options database which is initialized from command line arguments. Opendihu passes command line arguments on to PETSc such that this feature of PETSc can be used. 

Details can be found `here <https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetFromOptions.html>`_. For example by using the command line argument `-ksp_monitor` information about the PETSc solver are printed to stdout.
