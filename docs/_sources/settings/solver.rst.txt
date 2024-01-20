
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
        "dumpFilename": "",      # no filename means dump is disabled
        "dumpFormat": "default",
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


dumpFilename
~~~~~~~~~~~~~~
*Default: ""*

If this is set to a non-empty string, the system matrix and right hand side vector will be dumped before every linear solve. Dumping means writing to a file in a specific format, depending on
``dumpFormat``. 

Note, that there is also the possiblity of using :doc:`output_writer` s that output field variables. These also contain the right hand side vector, and any other variables of interest.

The value of ``dumpFilename`` is a filename prefix. The actual filename is of the form ``<dumpFilename>_matrix_#####.m``, where "matrix" is replaced by "rhs" for the right hand side vector and ##### is a 5-digit consecutive number, useful for timestepping schemes, where matrix and vector will be output for every time step. ".m" is the suffix, in this case for MATLAB files, otherwise ".txt" for ascii files.

A tip is to set it to ``"out/"``. This will create a directory called ``out``, the specified filename after this path is empty. But this is no problem, because the filename will only consist of the suffixes in this case.

dumpFormat
~~~~~~~~~~~~~~
*Default: "default"*

The format in which to export/dump data of matrices and vectors in the file. Possible values are:

- ``default``: Uses the default format provided by PETSc, usually human-readable, but not very nice. For sparse matrices this states the non-zero entries, for vectors it contains the values row-wise. This format is useful to copy-paste vectors to process them further in any other way.
- ``ascii``: Outputs dense vectors and matrices as ascii values, this is the best option to produce human-readable output. For a vector, the indices and values will be written row by row. For a matrix, the dense matrix including all zero entries will be written.
- ``matlab``: Generates files that can be used directly with MATLAB. Execute the ``*.m`` script files to load the variables in matlab. The variables will be named like the internal name in opendihu, e.g. for a Laplace problem in release mode ``rightHandSide`` and ``stiffnessMatrix``. The matrices are stored as sparse matrices in matlab.

The matlab option is also useful for higher-dimensional problems, whereas the other, ascii-based formats take a long time to write and should only be used with small problem sizes. The difference is that the matlab files contain the system matrix in a sparse format (i.e. only non-zero entries) and the ascii based files contain the whole matrix including all zeros.
 
