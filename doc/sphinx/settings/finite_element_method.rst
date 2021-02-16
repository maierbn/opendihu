FiniteElementMethod
====================

Description
------------

This class template is used to discretize the Laplace operator using the Finite Element Method. 
It assembles stiffness and mass matrices and can call a solver for the Laplace or Poisson problem.
Thus, it can be used directly to solve a Laplace or Poisson problem. It can also be used in a timestepping scheme to solve a diffusion equation.

C++ instantiation
-----------------

.. code-block:: c

  SpatialDiscretization::FiniteElementMethod<
    /*Mesh*/,
    /*BasisFunction*/,
    /*Quadrature*/,
    /*Equation*/
  >

Example:

.. code-block:: c
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<3>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  >
  
Mesh
^^^^^
The mesh to use for the discretization. This is one of 

* ``Mesh::StructuredRegularFixedOfDimension<D>``
* ``Mesh::StructuredDeformableOfDimension<D>``
* ``Mesh::UnstructuredDeformableOfDimension<D>``

where ``D`` is either ``1``, ``2`` or ``3``. See :doc:`/settings/mesh` for details.

BasisFunction
^^^^^^^^^^^^^^
The basis function/ansatz function/shape function to use for discretizing the field variables. This is one of 

* ``BasisFunction::LagrangeOfOrder<1>``
* ``BasisFunction::LagrangeOfOrder<2>``
* ``BasisFunction::Hermite``

where ``Hermite`` specifies cubic Hermite functions without scaling factors, i.e. :math:`C^0`, non-smooth discretizations.

Note that internally, the ``FunctionSpace::FunctionSpace`` class is build from ``Mesh`` and ``BasisFunction``.
  
Basis functions are defined in the directory `basis_function <https://github.com/maierbn/opendihu/tree/develop/core/src/basis_function>`_.

  
Quadrature
^^^^^^^^^^^
The quadrature scheme to use for computing the integral in the stiffness and mass matrices. For higher-dimensional quadrature, the tensor-product is taken. Possible values are

* ``Quadrature::Gauss<NumberGaussPoints>`` where ``NumberGaussPoints`` is the number of sampling points for the integration.
  
  Supported values are 1-8, 10, 12, 16, 20, 24, 64. 
  
  The Gauss quadrature scheme is capable of exactly integrating polynomials of degree :math:`2\cdot\text{NumberGaussPoints}-1`.
* ``Quadrature::ClenshawCurtis<NumberIntegrationPoints>`` for Clenshaw-Curtis integration. 
  
  It is capable of exactly integrating polynomials of degree :math:`\text{NumberIntegrationPoints}-1`. 
  
  Allowed values for ``NumberIntegrationPoints`` are 1-7 and 64.
* ``Quadrature::NewtonCotes<NumberIntegrationPoints>``.
 
  It is capable of exactly integrating polynomials of degree :math:`\lfloor(\text{NumberIntegrationPoints}+1)/2\rfloor \cdot 2-1`

  The following table lists the commonly known names of the schemes and their integration order.

  ==========  ==================  ================================
  name        integration points  exact integration of polynomials
  ==========  ==================  ================================
  rectangle   1                   1 
  trapezoid   2                   1 
  Simpson     3                   3
  3/8-rule    5                   5
  ==========  ==================  ================================
  
Quadrature schemes are defined in the directory `quadrature <https://github.com/maierbn/opendihu/tree/develop/core/src/quadrature>`_.
  
Equation
^^^^^^^^^^^
This specifies for which equation the mass and stiffness matrices should be computed. To directly use a ``FiniteElementMethod`` class with its ``run()`` method, the following values are possible:

* ``Equation::Static::Laplace`` for solving the Laplace equation:
  
  .. math::
    c\cdot Δ u = 0
    
* ``Equation::Static::GeneralizedLaplace`` for solving the Laplace equation:
  
  .. math::
    ∇\cdot(A∇u) = 0
    
  
* ``Equation::Static::Poisson`` for solving the Poisson equation:
  
  .. math::
    c\cdot Δ u = f
    
* ``Equation::Static::GeneralizedPoisson`` for solving the Poisson equation:
  
  .. math::
    ∇\cdot(A∇u) = f
    
  The right hand side :math:`f` can be given in the python settings.
    
For the following equations, the ``FiniteElementMethod`` class needs to be wrapped in other solver, e.g. time stepping schemes:

* ``Equation::Dynamic::IsotropicDiffusion`` for solving the diffusion equation,
  
  .. math::
    u_t - Δu = 0
    
* ``Equation::Dynamic::AnisotropicDiffusion`` for solving the diffusion equation,
  
  .. math::
    u_t - ∇\cdot A ∇u = 0
    
  
* ``Equation::Dynamic::DirectionalDiffusion`` for solving the diffusion equation,
  
  .. math::
    u_t - ∇\cdot A(v)∇u = 0,
    
  where the diffusion tensor, :math:`A`, depends on a direction field :math:`v`.
    

Equations are defined in the directory `equation <https://github.com/maierbn/opendihu/tree/develop/core/src/equation>`_. The prefactors :math:`c` and :math:`A` can vary over the domain.

  
Python Settings
---------------

The following keywords in the python dictionary are recognized:

.. code-block:: python

  "FiniteElementMethod" : {
    # <mesh>
    # <solver>
    # <output writer>
    
    "prefactor":          # type: double
    "rightHandSide":      # type: list of double
    "inputMeshIsGlobal":  # type: bool
    "diffusionTensor":    # type: list of double
    
    "nElements":          # type: integer 
    "physicalExtent":     # type: double
    "dirichletBoundaryConditions": # type: dict, {} 
    "neumannBoundaryConditions": # type: list, []
    "updatePrescribedValuesFromSolution": # type: bool
    "nodePositions":      # type: [[x,y,z], [x,y,z], ...]
    "elements":           # type: [[i1,i2,...], [i1,i2,...] ],
    "relativeTolerance":  # type: double
    "inputMeshIsGlobal":  # type: bool
    "slotName":           # type: string
    "OutputWriter":       # type: [{}, {}, ...]
  },

The items ``# <mesh>``, ``# <solver>`` and ``# <output writer>`` are placeholders for :doc:`/settings/mesh`, :doc:`/settings/solver` and :doc:`/settings/output_writer`.

.. _femesh:

<mesh>
^^^^^^^^^^^^^
The mesh defines the node positions and therefore the geometry. The number and location of the degrees of freedom (dofs) is dependent on the ``Mesh`` and ``BasisFunction`` classes in the C++ instantiation. Details on dofs and how to specify the geometry is given on the :doc:`/settings/mesh` page.

To specify properties of the mesh there are two possibilities:

1. Provide the ``meshName`` keyword and give the name of the mesh. A mesh with this name should be defined under the section ``Meshes`` of the config. This section is always top-level in the ``config`` dict, i.e. at ``config["Meshes"]``. For details, see :doc:`/settings/mesh`.

  An example for a 2D linear mesh is given below:
  
  .. code-block:: python

    config = {
      "Meshes": {
        "mesh0": {
          "nElements": [10, 10],
          "physicalExtent": [5.0,5.0]
          "inputMeshIsGlobal": True,    # This is always needed for the mesh
        }
      }
      "FiniteElementMethod" : {
        "meshName": "mesh0"
        
        # further options of FiniteElementMethod, e.g.:
        "dirichletBoundaryConditions": {0:0.0, -1:1.0},
        "prefactor": 1.0,
        "relativeTolerance": 1e-15,
        "rightHandSide": [0.0, 1.0, 2.0],
        "inputMeshIsGlobal": True   # This is needed at this point only, if a rightHandSide is given (for Poisson problem)
      }
    }
    
2. Directly specify the properties of the mesh under ``"FiniteElementMethod"``, for example:

  .. code-block:: python

    config = {
      "FiniteElementMethod" : {
        # options for the mesh
        "nElements": [10, 10],
        "physicalExtent": [5.0,5.0]
        "inputMeshIsGlobal": True,
        
        # further options of FiniteElementMethod, e.g.:
        "dirichletBoundaryConditions": {0:0.0, -1:1.0},
        "prefactor": 1.0,
        "relativeTolerance": 1e-15,
        "rightHandSide": [0.0, 1.0, 2.0],
      }
    }
    
  Note, that in this case, ``inputMeshIsGlobal`` applies to both the mesh data and the rightHandSide data, if specified
  
  
The first option is useful to reuse meshes that only need to be defined once. 

If a ``"meshName"`` is provided, all mesh properties specified under ``"FiniteElementMethod"`` other than ``"inputMeshIsGlobal"`` will be ignored.

<solver>
^^^^^^^^^^^^^
The solver is the solver of the linear system 

.. math::
  K u = f,

with stiffness matrix :math:`K`, vector of unknowns, :math:`u` and right hand side :math:`f`. This is needed when the Laplace or Poisson problem is solved (by calling ``run()`` of the object). Furthermore, in an explicit Euler timestepping of the diffusion equation we compute

.. math::
  u^{(t+1)} = u^{(t)} + dt M^{-1} K u_{t}

For this also the linear solver is used.

The specification of the solver can be given directly in-place or by specifying the name of a solver configuration that was given earlier, analoguous to <mesh>:

1. Provide the ``solverName`` keyword and give the name of the solver configuration. A solver configuration with this name should be defined under the section ``Solvers`` of the config. This section is always top-level in the ``config`` dict, i.e. at ``config["Solvers"]``. For details, see :doc:`/settings/solver`. 

  An example is given below:
  
  .. code-block:: python

    config = {
      "Solvers": {
        "linearSolver": {
        "relativeTolerance": 1e-15
        }
      }
      "FiniteElementMethod" : {
        "solverName": "linearSolver"
        
        # further options of FiniteElementMethod
      }
    }
  
2. Directly specify the solver properties under ``"FiniteElementMethod"``, for example:

  .. code-block:: python

    config = {
      "FiniteElementMethod" : {
        "relativeTolerance": 1e-15,
        
        # further options of FiniteElementMethod
      }
    }
  
The first option is useful when the same solver should be used for multiple classes.

If a ``"solverName"`` is provided, all solver properties specified under ``"FiniteElementMethod"`` will be ignored.

slotName
^^^^^^^^^^^^^^^^^

The `FiniteElementMethod` class exposes one slot that contains the solution variable. The option `slot name` specifies the name of the slot. Slot names are needed for coupling schemes to connect field variables between solvers. For details, see :doc:`output_connector_slots`.

<output writer>
^^^^^^^^^^^^^^^^^

This optionally specifies a list of output writers that can be used to output geometry field variables in various formats. For details, see :doc:`/settings/output_writer`

prefactor
^^^^^^^^^^
*Default: 1*

The prefactor is a scalar multiplier of the Laplace operator term, i.e. :math:`c` in :math:`c\cdot  Δu` or :math:`c\cdot∇ \cdot (A ∇u)`. 

It can be specified spatially varying, :math:`c(x)` by providing a list with as many entries as degrees of freedom (dofs). Depending on ``"inputMeshIsGlobal"``, you have to set as many entries as there are global or local dofs.

When using a composite mesh, you can also provide a list with as many entries as there are sub meshes. Then each value will be set in a sub mesh and :math:`c(x)` will be constant in the sub mesh.

rightHandSide
^^^^^^^^^^^^^^^
*Default: lists of zeros*

Here, the right hand side vector :math:`f` can be specified. 
Either as list of double values (e.g. ``[1.0, 0.0, 0.0, 2.0]``) or as dictionary (e.g. ``{0:1.0, 3:2.0}``). In the latter case all unspecified values are set to 0. 

In the dictionary case, negative indices are counted from the end. This means, e.g., that ``{-2:1.0, -1:1.0}`` sets the last two values to 1.0. This works for both values of ``inputMeshIsGlobal``, i.e. local and global indices (see below).

The given values correspond to global degrees of freedom, if ``"inputMeshIsGlobal": True`` is given (or omitted, because ``True`` is the default) or to local degrees of freedom of the process, if ``"inputMeshIsGlobal": False``.

dirichletBoundaryConditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See :doc:`boundary_conditions`.

neumannBoundaryConditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See :doc:`boundary_conditions`.

updatePrescribedValuesFromSolution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*Default:* ``False``

If this option is set to true, the values that are initially set in the solution field variable are used as the prescribed values at the dofs in `dirichletBoundaryConditions`.
The values that were given in `dirichletBoundaryConditions` have overridden by this. This is useful only if the `FiniteElementMethod` is part of a nested solver structure with a coupling and a timestepping scheme around it, where the solution value is updated in every iteration and the `solve()` gets called. Then the problem adjusts to update Dirichlet boundary conditions.o

inputMeshIsGlobal
^^^^^^^^^^^^^^^^^^
*Default:* ``True``

Together with ``rightHandSide`` it specifies whether the given values are interpreted as local values or global values in the context of a parallel execution on multiple processes. It has no effect for serial execution.

* If set to ``True``, values are interpreted as for serial execution. Then the same right hand side values should be given for all processes. Consequently, the program can be run with different numbers of processes with the same settings.
* If set to ``False``, the specified right hand side values are interpreted to be for the local portion of the own process. In parallel execution, each process has to get only its own range of values, which are typically different for each process. Only the non-ghost nodal values have to be given in the settings.

The advantage of the local specification is that each process only has to know its own portion of the whole problem. Internally there is no transfer of the local information to other processes. 
Thus, large problems can be computed with a high number of processes, where the global problem data would be too big to be stored by a single process.

To provide different values for different MPI ranks, the own MPI rank number can be retrieved in the python settings. The last two command line arguments that are available in the python settings script are the own MPI rank number and the total number of ranks.

The following example uses such a local specification of the right hand side. It sets the right hand side value of the last degree of freedom on the last MPI rank to 1.0 and all other values to 0.0.

.. code-block:: python

  # get own MPI rank number and number of MPI ranks
  rank_no = (int)(sys.argv[-2])
  n_ranks = (int)(sys.argv[-1])
  
  config = {
    "FiniteElementMethod" : {
      "inputMeshIsGlobal": False,
      "rightHandSide": {-1: 1.0} if rank_no == n_ranks-1 else {},
      
      # further options of FiniteElementMethod
      # ...
    }
  }

diffusionTensor
^^^^^^^^^^^^^^^^^^
For anisotropic diffusion, the diffusion tensor can be given as a list of double valus in row-major order. The diffusion tensor is always 2nd order (a square DxD matrix for dimension D). It is convenient to write the diffusion tensor with corresponding line breaks, as in the following example:

.. code-block:: python

  "FiniteElementMethod" : {
     "diffusionTensor": [
        8.93, 0, 0,
        0, 0.893, 0,
        0, 0, 0.893
      ], 
      
    # further options of FiniteElementMethod
    # ...
  }

The ``diffusionTensor`` can be specified spatially varying, :math:`A(x)` just like the ``prefactor``. This can be achieved by providing a list with as many entries as degrees of freedom (dofs). Depending on ``"inputMeshIsGlobal"``, you have to set as many items as there are global or local dofs. Note that one item is again a list of the entries of the matrix.

When using a composite mesh, you can also provide a list with as many items as there are sub meshes. Then each tensor will be set in a sub mesh and :math:`A(x)` will be constant in the sub mesh.

Properties
----------
* *Runnable*:   This class contains a ``run()`` method that solves the numerical problem. Therefore, this class can be used as the outermost solver of the instantiation in the ``main`` function.
* *Multipliable*:  This class can be wrapped by a ``MultipleInstances`` class.

Location in the code
--------------------
The source code files are in the directory `spatial_discretization/finite_element_method <https://github.com/maierbn/opendihu/tree/develop/core/src/spatial_discretization/finite_element_method>`_. 
The top-level class is in `finite_element_method.h <https://github.com/maierbn/opendihu/blob/develop/core/src/spatial_discretization/finite_element_method/finite_element_method.h>`_.


Examples
----------

The following example instantiates a finite element solver on a 3D structured regular grid, with linear Lagrange basis functions and Gauss quadrature with 3 Gauss points to solve the Laplace equation.

C++ main file

.. code-block:: c
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<3>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  >
  
python settings

.. code-block:: python

  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 5.0,
    "dirichletBoundaryConditions": {0:0.0, -1:1.0},
    "nodePositions": [[float(i)/n] for i in range(n+1)],
    "elements": [[i, i+1] for i in range(n)],
    "relativeTolerance": 1e-15,
    "inputMeshIsGlobal": True,
    "prefactor": 1.0,
    "OutputWriter" : []
  },
  

