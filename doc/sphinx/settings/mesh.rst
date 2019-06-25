
Mesh
=======

Meshes appear in the C++ code. A ``Mesh`` together with a ``BasisFunction`` define a ``FunctionSpace``. In the python settings, there are no function spaces, the *mesh* is everything that needs to be specified.

C++ class template Mesh
-------------------------

On the C++ side, there exist three mesh types, 

* ``Mesh::StructuredRegularFixedOfDimension<D>``,
* ``Mesh::StructuredDeformableOfDimension<D>``,
* ``Mesh::UnstructuredDeformableOfDimension<D>``,

with ``D`` equals to 1, 2 or 3.

:numref:`meshes` shows the different mesh types for linear and quadratic basis functions. As can be seen, the meshes all have quadrilateral elements, corresponding to lines in 1D, quadrilaterals in 2D and hexaeders in 3D.

.. _meshes:
.. figure:: meshes.svg
  :width: 50%

  Different Mesh types, here 2D meshes, for linear or Hermite basis functions (left column) and for quadratic Lagrange basis functions (right column)

The **Structured Regular Fixed** mesh is a structured rectilinear mesh with a fixed mesh width. The mesh width is the same in all dimensions. 
A StructuredRegularFixed mesh is fully defined by the number of elements in each coordinate direction and the mesh width.
The geometry of a mesh of this type cannot be changed during the computation, that means it is not possible to use this mesh for computing deformations.

One advantage of this mesh type is that the needed memory is independent of the number of elements and nodes, because only the mesh width and dimensions need to be stored. Therefore it is suited for large problem sizes.
A second advantage is that mass and stiffness matrices do not need to be integrated numerically if linear Lagrange basis functions are used. These matrices are instead assembled from precomputed stencils.

The **Structured Deformable** mesh is more general than the *Structured Regular Fixed* mesh. The node positions can have arbitrary locations.
However, the mesh is still structured, i.e. it has a fixed number of elements in the coordinate directions.

When computing elasticity problems, the geometry of the mesh (node positions) can change, i.e. the mesh can deform. 

This mesh type is the one most commonly used, when a given geometry should be discretized. For these two *structured* mesh types, the domain decomposition always uses planar cuts to subdivide the whole mesh into subdomains. Therefore, a subdomain again has always a fixed number of elements in each coordinate direction.

The **Unstructured** mesh is the most general mesh type. Contrary to the structured meshes, here the adjacency information can be defined arbitrarily and is not implicitely given by the mesh structured. 
The node positions need to be specified and can move during the computation, like with the *Structured Deformable* mesh.

The *Unstructured* mesh type can only be used with serial execution, i.e. no domain decomposition is implemented.

Node positions are always stored as points in :math:`\mathbb{R}^3`. Consequently, it is possible to define a 1D mesh embedded in the 3D space, for example for 1D muscle fibers in a 3D muscle geometry. Similarly, "bended" 2D meshes can be defined, like the 2D surface of a 3D muscle.
   
Python settings for Mesh
---------------------------
   
To specify a mesh in the python settings, depending on the mesh type, its size, node positions or adjacency infos have to be given.

In the python settings, all meshes should be defined at the beginning of the ``config`` dictionary under the ``"Meshes"`` key. 
The item ``"Meshes"`` is itself a dictionary, where the properties of every mesh are contained with a *mesh name* as key. The mesh name can be chosen abitrarily.
It is used to reference the mesh later, where it is needed, e.g. in a FiniteElementMethod object.

.. code-block:: python

  config = {
    "Meshes": {
      "mesh0": {
        "inputMeshIsGlobal": False,
        # further properties of mesh "mesh0" (see below for the description, depending on the type)
      },
      "arbitraryMeshName1": {
         # properties of this mesh
      }
    }
  }

StructuredRegularFixed
^^^^^^^^^^^^^^^^^^^^^^^ 
The **Structured Regular Fixed** mesh is completely specified by the number of elements in each coordinate direction and the physical extent.

.. code-block:: python

  "nElements": [nx, ny],      # example for a 2D mesh
  "physicalExtent": [2.5, 5.0],
  "inputMeshIsGlobal": True,

nElements
~~~~~~~~~~~~
*Default: D=1 (lines): 0, which means a degenerate element, D=2 or D=3: 1*

The number of elements of the mesh in the coordinate directions. For D=1, i.e. lines, it is a scalar non-negative integer value. For D=2 respective D=3 it is a list of 2 respective 3 non-negative integer values.
  
physicalExtent
~~~~~~~~~~~~~~~~
*Default: list of values 1.0*

The "size" of the mesh in physical units (e.g. meters if SI units are used), in the coordinate directions. This has to be a list of ``D``  positive double values.

Because the mesh width has to be constant in all coordinate directions, ``physicalExtent`` has to be a multiple of ``nElements``.

inputMeshIsGlobal
~~~~~~~~~~~~~~~~~~
*Default: ``True``*

Whether the values of ``nElements`` and ``physicalExtent`` describe the global domain (``True``) or the local subdomain (``False``) in parallel execution. See also the notes on :ref:`inputMeshIsGlobal` later.

StructuredDeformable
^^^^^^^^^^^^^^^^^^^^^^^ 
For specifying the **Structured Deformable** mesh. there are two possibilities: 

1. Specify ``nElements`` and ``physicalExtent``, like for a *StructuredRegularFixed* mesh. A rectilinear mesh is constructed, analogous to the *StructuredRegularFixed* mesh. 

  Note, that now the mesh widths does not need to be the same in every coordinate direction, so there is no restriction on the values of ``nElements`` and ``physicalExtent``.
  Again, the value of ``inputMeshIsGlobal`` applies.
  
2. Specify ``nElements`` and the node positions.

.. code-block:: python

  "nElements": [nx, ny],     # example for a 2D mesh
  "nodePositions": [[x,y,z], [x,y,z], ...], 
  "inputMeshIsGlobal": True,

nElements
~~~~~~~~~~~~
*Default: D=1 (lines): 0, which means a degenerate element, D=2 or D=3: 1*

The number of elements of the mesh in the coordinate directions. For D=1, i.e. lines, it is a scalar non-negative integer value. For D=2 respective D=3 it is a list of 2 respective 3 non-negative integer values.
  
nodePositions
~~~~~~~~~~~~~~~
Specify all node positions. There are two different formats:

1. A list of points where each point is a list with three entries ``[x,y,z]``. Even for lower dimensional meshes, ``D<3``, the node positions have three components.

  It is possible to define an embedded 1D or 2D manifold in the 3D domain. If this is not needed, the last entries can be set to 0 or omitted, as the default value for not specified components is 0.

2. The geometry values can also be given as consecutive array of [x,y,z,x,y,z,...] or [x,y,x,y,...] values.

  Then there is another property ``"nodeDimension"``, which is an integer from 1 to 3, with default value 3. This has to be set the number of dimensions that will be specified for each point in the consecutive array.

  Then, the node position values are provided in ``nodePositions`` as a list with ``nodeDimension`` double values per point, one point after each other.
  If ``nodeDimension`` is set to 1, ``nodePositions`` should be a list of the ``x`` values of the nodes, useful only for 1D meshes.
  If ``nodeDimension`` is set to 2, ``nodePositions`` should be a list with 2*number of nodes values, the x and y components of the node positions in consecutive order. Similar for ``nodeDimension=3``.

The order of the node positions proceeds through the entire structured mesh, with ``x`` advancing fastest, then the ``y`` index, then thet ``z`` index (if any). 
This means, e.g. for a 3D mesh, that starting from the first point at index :math:`(z,y,x)=(0,0,0)`, the next point is the one next to it in x-direction, i.e. :math:`(z,y,x)=(0,0,1)`,
then the next and so on until the line is full. Then the next line starts with :math:`(z,y,x)=(0,1,0)`, then :math:`(z,y,x)=(0,1,1)`, etc. 
After the x,y-plane is done, the next plane starts with :math:`(z,y,x)=(1,0,0)`.

For imagination see :numref:`coordinate_directions`.

.. _coordinate_directions:
.. figure:: coordinate_directions.svg
  :width: 80%

  Coordinate directions x,y,z and iterator/index names i,j,k for 2D and 3D meshes.
  
inputMeshIsGlobal
~~~~~~~~~~~~~~~~~~
*Default: ``True``*

Whether the values of ``nElements`` and the ``nodePositions`` describe the global domain (``True``) or the local subdomain (``False``) in parallel execution.

See also the notes on :ref:`inputMeshIsGlobal` later.

UnstructuredDeformable
^^^^^^^^^^^^^^^^^^^^^^^ 

For specifying an **Unstructured Deformable** mesh there are two options:

1. Using node positions and elements

2. Using *Exfiles*

.. _inputMeshIsGlobal:

inputMeshIsGlobal
^^^^^^^^^^^^^^^^^^^
It specifies whether the given values and degrees of freedom are interpreted as local values or global values in the context of a parallel execution on multiple processes. It has no effect for serial execution.
It applies to all values given as mesh properties, such as node positions, element and node numbers, the physicalExtent, the number of elements, etc.

* If set to ``True``, all specified values and degrees of freedom are interpreted with global indexing. In this case, the same values should be given on all processes. Consequently, the program can be run on different numbers of processes with the same settings.
* If set to ``False``, all specified values and degrees of freedom are interpreted to be for the local portion of the own process, only.
  In parallel execution, each process has to get only its own range of values, which are typically different on each process. 

  For example, the number of elements is only specified for the local portion. Opendihu will compute the global number of elements from the local numbers on all the processes.

To provide different values for different MPI ranks, the own MPI rank number can be retrieved in the python settings. 
The last two command line arguments that are available in the python settings script are the own MPI rank number and the total number of ranks.

The advantage of the local specification is that each process only has to know its own portion of the whole problem. Internally there is no transfer of the local information to other processes. 
Thus, large problems can be computed with a high number of processes, where the global problem data would be too big to be stored by a single process.

The following example uses such a local specification. It sets the right hand side value of the last degree of freedom on the last MPI rank to 1.0 and all other values to 0.0.

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
