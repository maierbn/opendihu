Boundary Conditions
===========================

Two types of boundary conditions are supported: Dirichlet and Neumann-type boundary conditions. They can be set in the settings for the following solvers:

* :doc:`FiniteElementMethod <finite_element_method>`
* :doc:`HyperelasticitySolver <hyperelasticity>`
* :doc:`DynamicHyperelasticitySolver <dynamic_hyperelasticity>`
* MuscleContractionSolver

Furthermore, in the :doc:`DynamicHyperelasticitySolver <dynamic_hyperelasticity>`, they can be changed over time.

Python settings
-----------------

In all of the above mentioned solvers, the syntax for specifying Dirichlet and Neumann-type boundary conditions is the same. The following keys are used:

.. code-block:: python

  config = {
    ...
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions":   neumann_bc,
    ...
  }
  
In order to not specify any boundary conditions, use an empty dict for Dirichlet boundary conditions and an empty list for Neumann boundary conditions:

.. code-block:: python

  "dirichletBoundaryConditions": {},
  "neumannBoundaryConditions": [],
    
Dirichlet Boundary Conditions
----------------------------------

Dirichlet boundary conditions fix several degrees of freedom to a prescribed value.
It has to be a Dict of ``{<dof no>: <value>}`` entries.

Negative dof nos are interpreted as counted from the end, i.e. -1 is the last dof, -2 is the second-last etc.
The ``<value>`` is a list of as many entries as the solution field has components. For example, for normal FiniteElementMethod, this is typically one component.
For the FiniteElementMethod that is used inside the 3D `QuasiStaticLinearElasticitySolver`, it is three components for the displacements in the three coordinate directions.
The same holds for `HyperelasticitySolver`. For the `DynamicHyperelasticitySolver`, the number of components is 6, comprising 3 displacements and 3 velocities.

``None`` can be used in place of a value to not prescribed a particular component.

The following example for a three-component FiniteElementMethod sets the degrees of freedoms (dof) with numbers `0, 5, 7` and the last degree of freedom to some values.
For dof 5, only the first component is prescribed.

.. code-block:: python

  dirichlet_bc = {
    0: [1,2],
    5: [3,None],
    -1: [5,6]
  }
  
  dirichlet_bc[7] = [7,7]
  
The example illustrates, how the dofs and values can be set. Either on construction of the dict within `{` and `}` or later with the `[]` operator.
  
Dirichlet boundary conditions are specified for dof numbers, not nodes, such that for Hermite it is possible to prescribe derivatives. For Lagrange ansatz functions, dof numbers are equivalent to node numbers.

The option `inputMeshIsGlobal` affects how the dirichlet bc configuration is interpret.
If `inputMeshIsGlobal` is set to ``True``, the numbering of the dofs is in *global natural* order. This means that the dofs are numbered fastest in `x`-direction then in `y`-direction, then in `z`-direction (Note for developers: this is different from the internal *global petsc* ordering of the actual memory layout in the local Petsc Vecs).

If `inputMeshIsGlobal` is set to ``False``, the specified dofs are interpreted as local numbers in the subdomain. Then you have to specify values **also for the ghost dofs**. This means that you have to specify prescribed nodal values for a node on every process whose subdomain is adjacent to that node.

The ghost dof numbers are after the non-ghost numbers. For example, consider the following mesh oft two linear elements, ``e1`` and ``e2`` on two ranks, ``r1`` and ``r2``:

.. code-block:: python

  dof numberings:
  local           global natural
  (e1)  (e2)      (e1)   (e2)
  1-3   2-3       4-5    6-7
  0-2   0-1       0-1    2-3
  r0     r1

Note how the left element has two ghost nodes, with local numbers 2 and 3 and how the local numbering is different from the right element which has no ghost nodes.

For **unstructured meshes**, the ordering of the dofs cannot be known at the time when the settings are parsed, because they depend on the mesh which could be read from ``*.ex`` files after the settings get parsed.
Therefore the ordering is special.
For every node there are as many values as dofs, in contiguous order.

Consider the following example for 2D Hermite, unstructured grid, 2x2 elements:

.. code-block:: python

  node numbering:
   6_7_8
  3|_4_|5
  0|_1_|2

  dof numbering:
   6_7_8
  2|_3_|5
  0|_1_|4

To specify du/dn = 0 at the left boundary in this example you would set:

.. code-block:: python
  
  bc[0*2+1] = 0, bc[3*2+1] = 0, bc[6*2+1] = 0

To specifiy u=0 on the bottom, you would set:

.. code-block:: python
  
  bc[0] = 0, bc[2] = 0, bc[4] = 0

For **composite meshes** the numbering proceeds through all sub mesh after each other. This means, numbers 0 to ``nDofsMesh0-1``, where ``nDofsMesh0`` is the number of dofs in the first submesh directly map to the dofs of the first submesh. Then the numbers ``nDofsMesh0`` to ``nDofsMesh0+nDofsMesh1-1`` map to the second sub mesh and so on. Note, that negative values therefore count from the end of the last submesh, i.e. ``-1`` specifies the last dof of the last submesh.

When the value to set is a vector, e.g. as mentioned for solid mechanics problems where displacements can be prescribed, specify a list of the components for each prescribed dof, e.g. ``[1.0, 2.0, 3.0]`` to set a Dirichlet boundary condition of :math:`\bar{u} = (1,2,3)^\top`. 
When not all components should be prescribed, replace the entry by ``None``, e.g. ``[None, 2.0, None]`` to only prescribe the y component. Or, if for the `DynamicHyperelasticitySolver` only displacements should fixed, used e.g. ``[0, 0, 0, None, None, None]``.


Neumann Boundary Conditions
----------------------------------

Neumann boundary conditions specify fluxes or traction forces for mechanics problems. These boundary conditions are specified for surfaces of elements. 

The given object is a list of dicts as shown below.

.. code-block:: python
  
  neumann_bc = [
    {
      "element": 0,               # local element no, negative values count from the end
      "face": "2-",               # face on which the neumann bc should act
      "constantVector": [1,0,0],  # specify only one of "constantVector", "constantValue" and "dofVectors"
      "constantValue": 0,
      "dofVectors":    {0:[1,0,0], 1:[2,0,0], 3:[2,1,0]},
    },
    {...}
  ]

* ``element`` is the local element number of the element which has the face for which to specify the Neumann boundary condition.
* ``face`` specifies the face on which the Neumann boundary condition will act. Possible values are "0-", "0+", "1-", "1+", "2-", "2+", where 0,1,2 stand for the x, y and z coordinate axis and "-" or "+" stand for negative or positive direction. Note that a 3D element has all 6 faces. A 2D element only has the faces "0-", "0+", "1-" and "1+". A 1D line element only has "0-" and "0+" which mean `left end point` and `right end point`.

The value to be prescribed can be given by either of three posibilities:

* ``constantVector``: This is a vector, useful e.g. for traction or forces. If for the :doc:`HyperelasticitySolver <hyperelasticity>` or :doc:`DynamicHyperelasticitySolver <dynamic_hyperelasticity>` the option `"divideNeumannBoundaryConditionValuesByTotalArea"` is set to `True`, this vector is interpreted as a total force and will be scaled down automatically to reflect the actual surface size. If this option is `False`, the vector specifies a constant surface traction of the whole surface, as one would expect. The direction is always given in reference configuration. Note that you are responsible for the direction of the vector, if it points inwards or outwards of the domain.
* ``constantValue``: For problems with only 1 component, this is the natural way to specify an outward flux. For problems with more components, i.e. mechanics, the constant value will be a factor to the outward normal. By this it is easily possible to specify loads orthogonal to the surface.
* ``dofVectors``: This gives the most flexibility in specifying the values. It is a dict where the key is the node/dof number on the *surface element* and the value specifies the traction vector on that node. In the given example, only the nodes 0, 1 and 3 have a Neumann boundary condition value applied.

Especially for the mechanics problems, it is often convenient to use Python list comprehension to construct the ``neumann_bc`` object:

.. code-block:: python
  
  nx = ... # number of elements in x direction
  ny = ... # number of elements in x direction
  nz = ... # number of elements in x direction
  k = nz-1   # top element
  
  # the following specifies a constant surface load of 1 pointing upwards and acting on the top surface of the whole 3D box
  neumann_bc = [{"element": k*nx*ny + j*nx + i, "constantVector": [1,0,0], "face": "2+"} for j in range(ny) for i in range(nx)]

