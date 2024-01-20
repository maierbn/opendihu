Hyperelasticity
====================

Description
------------

This class is a nonlinear solver for hyperelastic solid mechanics problems in 3D using structured meshes. Incompressible materials are implemented by mixed-formulation Finite Elements with Taylor-Hood elements, which use quadratic shape functions for displacements (and velocities for dynamic problems) and linear shape functions for the pressure.
For compressible materials, the pressure is computed from a constitutive equation and does not need an extra discretization.

This solver computes the static problem, for the dynamic problem, use the :doc:`dynamic_hyperelasticity` solver.


C++ instantiation
-----------------

.. code-block:: c

  SpatialDiscretization::HyperelasticitySolver<
    Material
  >
  // or:
  SpatialDiscretization::HyperelasticitySolver<Material, true>   // default, same as without "false"
  SpatialDiscretization::HyperelasticitySolver<Material, false>

If the second tempate parameter is ``true``, additionally the PK1 stress :math:`P` (instead of only the PK2 stress :math:`S`) and the deformation gradient :math:`F` will be contained in the output files. However, since :math:`P` is unsymmetric it will be more data (9 values per dof). If the second parameter is set to ``false``, the output files will be significantly smaller. 
The first template parameter ``Material`` is a class that describes the used constitutive equations at compile time.

Specification of the Material
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following classes are pre-defined for the material:

.. code-block:: c

  Equation::SolidMechanics::MooneyRivlinIncompressible3D
  Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressible3D
  Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D
  Equation::SolidMechanics::HyperelasticTendon
  
You can look at `core/src/equation/mooney_rivlin_incompressible.h <https://github.com/maierbn/opendihu/blob/develop/core/src/equation/mooney_rivlin_incompressible.h>`_ to see how they are defined.

Alternatively, a new formulation can directly be defined in the main cpp file. For this, the ``struct Material`` that is used as template parameter in the ``HyperelasticitySolver`` class has to be defined.
It contains three parts: 

* Some settings by defining ``bool`` constants
* Definition of material parameters. Values for these are set in the python settings.
* The strain energy function for the material.

The following is an example for this.

.. code-block:: c

  #include <Python.h>
  #include <iostream>
  #include <cstdlib>

  #include <iostream>
  #include "easylogging++.h"

  #include "opendihu.h"

  // define material
  struct Material : Equation::SolidMechanics::HyperelasticityBase
  {
    static constexpr bool isIncompressible = false;      //< if the formulation is incompressible, then, strainEnergyDensityFunctionVolumetric will not be considered
    static constexpr bool usesFiberDirection = false;    //< if the decoupled form uses the 4th or 5th invariants, Ibar4, Ibar2, this means it is an anisotropic material
    static constexpr bool usesActiveStress = false;      //< if the value of an active stress term will be added to the stress
    
    // material parameters
    static constexpr auto c1 = PARAM(0);     //< material parameter
    static constexpr auto c2 = PARAM(1);     //< material parameter
    static constexpr auto kappa = PARAM(2);  //< material parameter, bulk modulus

    static constexpr int nMaterialParameters = 3;  //< number of material parameters

    //! the isochoric part of the decoupled strain energy function, Ψ_iso(Ibar1,Ibar2,Ibar4,Ibar5), in terms of the reduced invariants
    static const auto constexpr strainEnergyDensityFunctionIsochoric = c1*(Ibar1 - INT(3)) + c2*(Ibar2 - INT(3));

    //! the volumetric part of the decoupled strain energy function, Ψ_vol(J), only used for compressible formulation (isIncompressible == false)
    static constexpr auto G = INT(1)/INT(4) * (pow(J,INT(2)) - INT(1) - INT(2)*ln(J));        // Holzapfel p.245
    static const auto constexpr strainEnergyDensityFunctionVolumetric = kappa * G;

    //! coupled form of the strain energy function, Ψ(I1,I2,I3), as alternative to the two decoupled functions
    static const auto constexpr strainEnergyDensityFunctionCoupled = INT(0);

    //! another coupled form of the strain energy function, Ψ(C), dependent on right Cauchy Green tensor, C
    //! it must only depend on variables C11, C12, C13, C22, C23, C33.
    static const auto constexpr strainEnergyDensityFunctionCoupledDependentOnC = INT(0);
  };

  int main(int argc, char *argv[])
  {
    // Solves nonlinear hyperelasticity (Mooney-Rivlin) using the built-in solver
    
    // initialize everything, handle arguments and parse settings from input file
    DihuContext settings(argc, argv);
    
    
    // define problem
    SpatialDiscretization::HyperelasticitySolver<Material> problem(settings);
    
    // run problem
    problem.run();
    
    return EXIT_SUCCESS;
  }
  
Every new material has to inherit from ``Equation::SolidMechanics::HyperelasticityBase``, which defines the following symbols:

.. _baseclass:

.. code-block:: c

  // reduced invariants, arguments of `strainEnergyDensityFunctionIsochoric`
  static constexpr auto Ibar1; //< 1st reduced or modified strain invariant Ibar1 = tr(Cbar) = J^{-2/3}*I_1
  static constexpr auto Ibar2; //< 2nd reduced or modified strain invariant Ibar2 = 1/2 (tr(Cbar)^2 - tr(Cbar^2)) = J^{-2/3}*I_2
                                        // Note, there is no 3rd reduced or modified strain invariant needed, Ibar3 = det(Cbar) = 1 (incompressibility)
  static constexpr auto Ibar4; //< 4th reduced or modified strain invariant Ibar4 = a0•C a0
  static constexpr auto Ibar5; //< 5th reduced or modified strain invariant Ibar4 = a0•C^2 a0
  static constexpr auto lambda = sqrt(Ibar4);    //< fiber stretch, helper variable that can also be used in `strainEnergyDensityFunctionIsochoric`

  // volume factor, argument of `strainEnergyDensityFunctionVolumetric` (only for compressible material)
  static constexpr auto J;     //< volume factor, J = det(F), only for compressible material (otherwise it is 1)

  // invariants, arguments of `strainEnergyDensityFunctionCoupled`
  static constexpr auto I1;    //< 1st strain invariant I1 = tr(C)
  static constexpr auto I2;    //< 2nd strain invariant I2 = 1/2 (tr(C)^2 - tr(C^2))
  static constexpr auto I3;    //< 3rd strain invariant I3 = det(C) = J^2

  // components of the right Cauchy Green tensor, arguments of `strainEnergyDensityFunctionCoupledDependentOnC
  static constexpr auto C11;    //< entry C11 of the right Cauchy Green tensor, C
  static constexpr auto C12;    //< entry C12 = C21 of the right Cauchy Green tensor, C
  static constexpr auto C13;    //< entry C13 = C31 of the right Cauchy Green tensor, C
  static constexpr auto C22;    //< entry C22 of the right Cauchy Green tensor, C
  static constexpr auto C23;    //< entry C23 = C32 of the right Cauchy Green tensor, C
  static constexpr auto C33;    //< entry C33 of the right Cauchy Green tensor, C

  static constexpr auto a1;     //< entry a0_1 of the fiber direction, a0
  static constexpr auto a2;     //< entry a0_2 of the fiber direction, a0
  static constexpr auto a3;     //< entry a0_3 of the fiber direction, a0

  static constexpr auto I4;     //< non-reduced 4th strain invariant, I4 = a0•C a0
  
These symbols are to be used as the parameters to the strain energy functions and are, thus, available in the material description class.
Only some particular symbols can be used in some terms of the strain energy function.

In the following, the three parts of a custom material are explained.

Specification of the options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Currently, three options have to be defined.

* ``isIncompressible`` specifies if the material is incompressible. This implies ``J=1`` and ``I3=1``. 
* ``usesFiberDirection``: if this is set, the material can be anisotropic and the invariants ``Ibar4`` and ``Ibar5`` are available. Furthermore, a fiber direction will be determined from the given fiber meshes in the python settings.
* ``usesActiveStress``: if this is true, the value of thet active stress will be added to the normal stress

Specification of the parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Parameters are constants that can appear in the formulation of the strain energy function. Their value has to be set in the python settings.

The parameter section looks like the following.

.. code-block:: c

    // material parameters
    static constexpr auto c1 = PARAM(0);     //< material parameter
    static constexpr auto c2 = PARAM(1);     //< material parameter
    static constexpr auto kappa = PARAM(2);  //< material parameter, bulk modulus

    static constexpr int nMaterialParameters = 3;  //< number of material parameters

Any number of parameters can be specified and the names are custom. (The specifiers already used in the ``HyperelasticityBase`` class can, of course, not be used).
The parameters are assigned the macro ``PARAM(i)`` where ``i`` is a consecutively increasing number from 0.
The number of parameters in ``nMaterialParameters`` has to be correct. This is the number of values that are expected in the python settings ``materialParameters``.
The order of the values in the python settings is given by the ``PARAM`` macros.

.. _strain_energy_function:

Specification of the strain energy function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The strain energy function, ψ, is the constitutive equation that connects deformation with stress response. For the 2nd Piola-Kirchhoff, :math:`S`, and the right Cauchy-Green tensor, :math:`C=F^\top\,F`, the following holds:

.. math::
  S = 2\dfrac{∂ψ}{∂C}

For a hyperelastic material, the strain energy can always be formulated in terms of invariants. The implemented functionality comprises 5 invariants. The first three specifiy isotropic material behaviour and the 4th and the 5th depend on a fiber direction.

.. math::
  
  I_1 &= tr(C),\\
  I_2 &= \dfrac12 \big(tr(C)^2 - tr(C^2)\big),\\
  I_3 &= det(C) = J^2,\\
  I_4 &= a_0 \cdot C a_0,\\
  I_5 &= a_0 \cdot C^2 a_0,\\

where :math:`C` is the right Cauchy-Green tensor, :math:`J = det F` is the volume ratio or the determinant of the deformation gradient, :math:`F = grad_X χ` and :math:`a_0` is the fiber direction.

Often, the strain energy is given in a decoupled formulation. It is formulated in terms of quantities that are split into volumetric and isochor (with constant volume) contributions.  We have the volume-preserving measures

.. math::

  \bar{F} = J^{-1/3}\,F, \quad \bar{C} = J^{-2/3}\,C

Then, we use the modified or reduced invariants

.. math::
  
  \bar{I}_1 &= tr(\bar{C}) &= J^{-2/3}\,I_1,\\
  \bar{I}_2 &= \dfrac12 \big(tr(\bar{C})^2 - tr(\bar{C}^2)\big) &= J^{-4/3}\,I_2,\\
  \bar{I}_3 &= det(\bar{C}) = 1.

The general form in which the strain energy function can be specified consists of the following 4 summands.

.. math::
  Ψ = Ψ_{iso}(\bar{I}_1, \bar{I}_2, \bar{I}_4, \bar{I}_5) + Ψ_{vol}(J) + Ψ(I_1,I_2,I_3) + Ψ(C,a_0,I4)
  
Every summand can be set to constant 0 if not needed (``INT(0)`` in the C++ code).

In order to use a decoupled formulation, specify :math:`Ψ_{iso}(\bar{I}_1, \bar{I}_2, \bar{I}_4, \bar{I}_5)` and :math:`Ψ_{vol}(J)` for compressible materials or only :math:`Ψ_{iso}(\bar{I}_1, \bar{I}_2, \bar{I}_4, \bar{I}_5)` for incompressible materials.

To use a coupled formulation, use :math:`Ψ(I_1,I_2,I_3)`. Though the strain energy function can always be formulated in terms of the invariants, some literature only provides a formulation in terms of the right Cauchy-Green tensor, :math:`C` and the fiber direction, :math:`a_0`. In this case, the function :math:`Ψ(C,a_0,I4)` can be specified. For the last function, note that :math:`I_4` is available as an abbreviation for :math:`a_0 \cdot C a_0`.

The available summands of :math:`Ψ` also depends on the options that were set in the first part of the material structure. For incompressible material, i.e. if ``isIncompressible == true``, we have the following form:

.. math::
  Ψ = Ψ_{iso}(\bar{I}_1, \bar{I}_2, \bar{I}_4, \bar{I}_5) + Ψ(I_1,I_2,I_3) + Ψ(C,a_0)
  
If ``usesFiberDirection == false`` there are no 4th and 5th invariants:

.. math::
  Ψ = Ψ_{iso}(\bar{I}_1, \bar{I}_2) + Ψ_{vol}(J)  + Ψ(I_1,I_2,I_3) + Ψ(C,a_0,I4)
  
The 4 functions :math:`Ψ_{iso}(\bar{I}_1, \bar{I}_2, \bar{I}_4, \bar{I}_5)` :math:`Ψ_{vol}(J)`, :math:`Ψ(I_1,I_2,I_3)` and :math:`Ψ(C,a_0)` are given by the following 4 symbols that need to be defined in the material struct:

.. code-block:: c

    static const auto constexpr strainEnergyDensityFunctionIsochoric = INT(0);      // parameters: Ibar1,Ibar2,Ibar4,Ibar5,lambda (=sqrt(Ibar4))
    static const auto constexpr strainEnergyDensityFunctionVolumetric = INT(0);     // parameters: J
    static const auto constexpr strainEnergyDensityFunctionCoupled = INT(0);        // parameters: I1,I2,I3
    static const auto constexpr strainEnergyDensityFunctionCoupledDependentOnC = INT(0);  // parameters: C11, C12, C13, C22, C23, C33, a1, a2, a3, I4
  
The equations need to be specified according to the syntax of the `SEMT library <https://github.com/st-gille/semt>`_. 
Normal operators such as ``+``, ``*``, ``sqrt``, ``ln`` and ``pow`` can be used to combine the parameters given under :ref:`the base class<baseclass>`. 
Whenever an integer constant needs to be used, wrap it in ``INT()``, e.g. ``INT(5)``. Other factors that are no whole numbers cannot be used directly. They have to be defined as material parameter and their value is then set in the python settings.

It is also possible to define helper functions that are reused later. This can be done with the type ``static constexpr auto``.

An example for the incompressible Mooney-Rivlin material is given below:

.. code-block:: c
  
  static const auto constexpr strainEnergyDensityFunctionIsochoric
    = c1*(Ibar1 - INT(3)) + c2*(Ibar2 - INT(3));
  
An example for an incompressible material that uses a helper function is given here:

.. code-block:: c
  
  static constexpr auto d = INT(2)*(c1 + INT(2)*c2);
  
  static const auto constexpr strainEnergyDensityFunctionCoupled 
    = c*pow(sqrt(I3) - INT(1), INT(2)) - d*ln(sqrt(I3)) + c1*(I1 - INT(3)) + c2*(I2 - INT(3));

  
Python settings
-----------------

The following shows all possible options. The meaning can be learned from the comments.

.. code-block:: python

  "HyperelasticitySolver": {
    "durationLogKey":             "duration_mechanics",         # key to find duration of this solver in the log file
    
    "materialParameters":         material_parameters,          # material parameters of the Mooney-Rivlin material
    "displacementsScalingFactor": 1.0,                          # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
    "residualNormLogFilename":    "log_residual_norm.txt",      # log file where residual norm values of the nonlinear solver will be written
    "slotNames":                  ["ux", "uy", "uz"],           # (optional) slot names of the data connector slots, there are three slots, namely the displacement components ux, uy, uz
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
    "snesMaxIterations":          100,                           # maximum number of iterations in the nonlinear solver
    "snesRelativeTolerance":      1e-5,                         # relative tolerance of the nonlinear solver
    "snesLineSearchType":         "l2",                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesAbsoluteTolerance":      1e-5,                         # absolute tolerance of the nonlinear solver
    "snesRebuildJacobianFrequency": 1,                          # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
    
    #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),          # dump system matrix and right hand side after every solve
    "dumpFilename":               "",                           # dump disabled
    "dumpFormat":                 "default",                     # default, ascii, matlab
    
    #"loadFactors":                [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
    #"loadFactors":                [0.5, 1.0],                   # load factors for every timestep
    #"loadFactors":                list(np.logspace(-3,0,4)),    # load factors, equally spaced in log space: (1e-3, 1e-2, 1e-1, 1)
    "loadFactors":                [],                           # no load factors, solve problem directly
    "loadFactorGiveUpThreshold":  4e-2,                         # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the progression between two subsequent load factors gets smaller than this value, the solution is aborted.
    "scaleInitialGuess":          False,                        # when load stepping is used, scale initial guess between load steps a and b by sqrt(a*b)/a. This potentially reduces the number of iterations per load step (but not always).
    "nNonlinearSolveCalls":       1,                            # how often the nonlinear solve should be called
    
    # boundary and initial conditions
    "dirichletBoundaryConditions": elasticity_dirichlet_bc,             # the initial Dirichlet boundary conditions that define values for displacements u
    "neumannBoundaryConditions":   elasticity_neumann_bc,               # Neumann boundary conditions that define traction forces on surfaces of elements
    "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
     
    "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
    "constantBodyForce":           constant_body_force,                 # a constant force that acts on the whole body, e.g. for gravity
    
    "dirichletOutputFilename":     "out/"+scenario_name+"/dirichlet_boundary_conditions_tendon",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    # define which file formats should be written
    # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
    "OutputWriter" : [
      
      # Paraview files
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      
      # Python files and callback
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/all/"+scenario_name, "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonCallback", "outputInterval": 1, "filename": "out/all/"+scenario_name, "callback": handle_result_hyperelasticity, "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    ],
    # 2. additional output writer that writes also the hydrostatic pressure
    "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
      "OutputWriter" : [
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
    # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  },

durationLogKey
^^^^^^^^^^^^^^^^
A key under which the duration for this solver is stored in the log file.

`materialParameters`
^^^^^^^^^^^^^^^^^^^^^^^
A list of material parameters, must match the number of parameters in the material.

displacementsScalingFactor"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A scaling factor for the displacements that will be written to the output files. This is mainly for debugging.
Only set this to something other than 1 to increase the visual appearance for very small displacements.

residualNormLogFilename
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A txt log file where the residual norm values of the nonlinear solver will be written to. 

The progression of the residual norm over number of iterations can be visualized using ``plot_residual_norm.py``.

slotNames
^^^^^^^^^^^^
(optional) A list of names for the data connector slots. The slot names are used for connecting the slots to other solvers, i.e., when the displacement results should be reused by another solver. 
Each slot name should have <= 10 characters. See :doc:`output_connector_slots` for more details.

The key `slotNames` can also be omitted if the slots should not be reused. Note that these slotNames are not needed if the Hyperelasticity solver is contained in a :doc:`muscle_contraction_solver`.

`useAnalyticJacobian` and `useNumericJacobian`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Whether to use the analytically computed jacobian matrix in the nonlinear solver (fast) or the numerically computed jacobian matrix in the nonlinear solver (slow). This only works with non-nested matrices, if both numeric and analytic are enabled, it uses the analytic for the preconditioner and the numeric as normal jacobian.
  
dumpDenseMatlabVariables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow). This is mainly for debugging.
If `useAnalyticJacobian`, `useNumericJacobian` and `dumpDenseMatlabVariables` are all three set to ``True``, the analytic and numeric Jacobian matrices will get compared to see if there are programming errors for the analytic jacobian. Use this only for very small problems (like 5 elements)

meshName
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The mesh to use, this mesh has to use quadratic Lagrange basis functions. See :doc:`mesh` how to specify meshes.

inputMeshIsGlobal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This refers to the specification of the boundary conditions. Indicates whether the numberings used in the BCs is interpreted as global or local numbers. Note, that the mesh can be specified independently, i.e., it is possible to have the mesh specification in local numberings and the boundary conditions in global numberings.

fiberMeshNames
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Fiber meshes that will be used to determine the fiber direction, used for anisotropic materials

fiberDirection, fiberDirectionInElement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If fiberMeshNames is empty, directly set the constant fiber direction. Either use `fiberDirection` to specify a direction in the global coordinate system, or use `fiberDirectionInElement` to specify the direction in the element coordinate system. The direction should be a vector, e.g., ``[0,0,1]``

Nonlinear Solver
^^^^^^^^^^^^^^^^^^^^^^^
The following parameters can be given to specify the nonlinear solver:

.. code-block:: python

  "relativeTolerance":          1e-5,                         # 1e-10 relative tolerance of the linear solver
  "absoluteTolerance":          1e-10,                        # 1e-10 absolute tolerance of the residual of the linear solver       
  "solverType":                 "preonly",                    # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
  "preconditionerType":         "lu",                         # type of the preconditioner
  "maxIterations":              1e4,                          # maximum number of iterations in the linear solver
  "snesMaxFunctionEvaluations": 1e8,                          # maximum number of function iterations
  "snesMaxIterations":          100,                           # maximum number of iterations in the nonlinear solver
  "snesRelativeTolerance":      1e-5,                         # relative tolerance of the nonlinear solver
  "snesLineSearchType":         "l2",                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
  "snesAbsoluteTolerance":      1e-5,                         # absolute tolerance of the nonlinear solver
  "snesRebuildJacobianFrequency": 1,                          # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
  
  "dumpFilename":               "",                           # dump disabled 
  "dumpFormat":                 "default",                    # default, ascii, matlab

Details, e.g., about `dumpFilename` can also be found under :doc:`solver`.

Note that the top-level ``regularization`` option also influences the error of the solution, see next section.

regularization
^^^^^^^^^^^^^^^^^^
When using non-cartesian meshes, some elements may have interior angles of almost 180 degrees. In such elements, the inversion of the mapping from element coordinates to world coordinates is nearly singular.
Inverting the jacobian matrix of this mapping is poorly conditioned and gives large errors. As a remedy, the following regularization can be specified. If the absolute determinant of the matrix is below a given threshold `tol`, a constant `eps*h` is added on the diagonal of the matrix. Thus, the problem becomes better conditioned. 
The constant depends on the mesh width `h`. For :math:`h \to 0`, the regularization vanishes and the formulation converges to the correct solution.

This regularization is applied to all matrix inversions (not only for the described mapping between element and world coordinates). The parameters `tol` and `eps` can be specified by the ``"regularization"`` option as a list ``[tol, eps]``.
The option can be set to ``None``, to disable the regularization (equivalent to ``[0,0]``). If not specified, the default regularization is `[1e-2,1e-1]`.

.. code-block:: python

  config = {
    "regularization":           [1e-5,1e-3],                  # None to disable or [tol,eps]. if the absolute determinant |J|=|det(F)| within an element is below tol, add eps*I to the matrix, to regularize the inversion of the nearly singular matrix
    #"regularization":           None,                         # None to disable or [tol,eps]. if the absolute determinant |J|=|det(F)| within an element is below tol, add eps*I to the matrix, to regularize the inversion of the nearly singular matrix
  }

loadFactors
^^^^^^^^^^^^^^^^^^^^^^^

The load factors to solve the static problem. This should be a list of factors between 0 and 1 with the last factor 1. The loads will be scaled with these factors. After a solution with a factor was solved, the next solution uses the previous solution as initial values. Thus, it is possible to solve badly conditioned problems by increasing the load step by step.

Examples for load factors:

.. code-block:: python
  
  [0.1, 0.2, 0.35, 0.5, 1.0],
  list(np.logspace(-2,0,10)),   # use 10 equidistant load factors in [0,1] in log space
  [0.5, 1.0],
  [],                           # no load factors, solve problem directly

loadFactorGiveUpThreshold
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When a computation diverges, it is automatically retried with a load factor that is half as large as the failed load step. 
This can lead to the load factors getting smaller and smaller without any successfull solution. In such a case it is desirable to abort the computation.
If the progression between two subsequent load factors gets smaller than this threshold value, the solution is finally considered diverged and the computation continues with the next solver.
    
scaleInitialGuess
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(default: False) After a load step has been computed, scale the resulting solution that is used as the initial guess for the next load step.
If the previous load factor is :math:`b` and the next load factor is :math:`b`, the naive way would be to scale by the factor :math:`b/a` (divide by old factor, multiply by new factor).
However, this would lead to an overshoot, as the material is approximated linearly but the real model is nonlinear. 
Instead, we scale scale by :math:`\sqrt{ab}/a`. This corresponds to the geometric mean between the old load factor and the new load factor.

This scaling usually reduces the initial residual. Nevertheless, the number of iterations is sometimes higher, maybe because the prediction led to a worse area in the definition space of the model.
  
Note, this option is different from `extrapolateInitialGuess`, which only applies to dynamic problems and uses information from the last timestep. The option `scaleInitialGuess` uses information from the previous load step and is indepent of whether the problem is static or dynamic.
  
nNonlinearSolveCalls
^^^^^^^^^^^^^^^^^^^^^^^

How often the same static problem should be solved. This should be set to 1, because it makes no sense to solve the same problem multiple times. It originates from the Chaste documentation, where they observed different solutions after the first solve (which doesn't make sense).


Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^
Boundary conditions are specified with the keys ``dirichletBoundaryConditions``, ``neumannBoundaryConditions`` and ``divideNeumannBoundaryConditionValuesByTotalArea``.
Refer to :doc:`boundary_conditions` how to specify boundary conditions.

``divideNeumannBoundaryConditionValuesByTotalArea`` specifies if the given Neumann boundary condition values under ``neumannBoundaryConditions`` are total forces or surface loads. If ``True`` the values are surface loads and will be scaled by the surface area of all elements where Neumann BC are applied. The unit is then `N/cm^2`. If ``False``, the values are treated as normal Neumann boundary condition values, i.e. nodal force values with unit `N`.

dirichletOutputFilename
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable. This is for debugging the Dirichlet boundary condition nodes

Initial values
^^^^^^^^^^^^^^^^^^^
The initial values are given by ``initialValuesDisplacements`` and ``initialValuesVelocities``. A list of entries for all dofs is required, as vector of values for every node: ``[[node1-x,y,z], [node2-x,y,z], ...]``

extrapolateInitialGuess
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities (from the previous timestep). 
This is faster and should be set to ``True``.

constantBodyForce
^^^^^^^^^^^^^^^^^^^^^^^^^^^
A constant force that acts on the whole body, e.g. for gravity. The units are ``cm/ms^2`` It should be a 3d vector (list of 3 entries):

.. code-block:: python
  
  constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
  
OutputWriters
^^^^^^^^^^^^^^^^^^
There are different types of output writers that output different variables.

* ``"OutputWriter"``: 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
* ``"pressure"``: 2. additional output writer that writes also the hydrostatic pressure
* ``"LoadIncrements"``: 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written


