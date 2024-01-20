
Features
============

This page lists the functionality that is currently implemented. More detailed information is given on the :doc:`/settings` pages.

Equations
-----------
OpenDiHu currently implements the following model equations:

* Laplace and Poisson
* Diffusion, istropoic and anisotropic
* Monodomain: Reaction-diffusion, where the reaction term is given in `CellML <https://www.cellml.org/>`_ description
* Bidomain equation
* Multidomain equation
* Multi-physics combinations of the above
* Functionality to estimate and trace fibers from streamlines of a laplacian flow through a muscle volume
* Static linear elasticity
* Hyperelasticity, static and dynamic formulations
* Pre-defined hyperelastic materials such as Mooney-Rivlin or user-defined strain energy functions

Boundary and initial conditions
-------------------------------------

* Dirichlet and Neumann-type boundary conditions are supported.
* For transient problems, initial conditions have to be specified.

Discretization and Mesh types
------------------------------------
* For the time domain, there are explicit and implicit Euler and Heun method implemented, Crank-Nicolson, as well as a Godunov and Strang splitting scheme. 

  To learn more about the time stepping schemes, see :doc:`/settings/timestepping_schemes_ode`, the splitting schemes are explained at :doc:`/settings/splitting`.
* There are three implemented mesh types which each have variants for 1D, 2D and 3D.

  * **StructuredRegularFixedOfDimension<D>**: 
  
    This is a rectilinear grid with fixed mesh width in all dimensions. To store a mesh, only the number of elements in the `D` dimensions and the mesh width or extent is needed. 
    
    This is computationally fast and memory efficient. The stiffness and mass matrices for the laplace operator are assembled by precomputed stencils. This mesh type is mainly used for validation of the other mesh types and for examples that have an analytical solution.
  * **StructuredDeformableOfDimension<D>**: 
  
    This structured mesh can have arbitrary node positions, with a fixed number of elements in each dimension. 
    
    Internally, there is no need to store extra adjacency information, only the node positions and number of elements in the dimensions is stored. Therefore, using this mesh type is computationally efficient. This type of mesh can deform over time, it is suitable for structured mechanics calculations.

  * **UnstructuredDeformableOfDimension<D>**: 
    
    This is the most general mesh type. The elements can arbitrarily share nodes. Supported are also *versions* like in OpenCMISS where there can be multiple degrees of freedom at a node, to model discontinuities. Such meshes can be natively constructed from `exfiles` which is a file format defined in the OpenCMISS ecosystem.

* The following basis functions are available: linear and quadratic Lagrange and cubic Hermite. Unit tests cover all combinations of mesh types in 1D, 2D, 3D with the basis function types, in serial and parallel execution.

More details on the meshes are given on the :doc:`/settings/mesh` page.

Input/Output
---------------

* The OpenDiHu framework has a Python 3 interpreter included and defines its own python dictionary-based config file format. Input for meshes can be specified by python lists containing node positions, element-node mappings, number of elements, etc., depending on the mesh type. This gives the user flexibility as to how to import and store custom mesh information.
* As a special case, `exnode` and `exelem` file types can be parsed to construct a mesh of type `UnstructuredDeformableOfDimension<D>`.
* OpenDiHu supports 5 output formats for simulation results:

  * **Python file:** 
  
    The data is exported as python dictionary. The output file can be configured to be human readable (ascii), or binary (using the python `pickle` package). 
    
    This type of output is the most memory efficient and easy to handle output format. OpenDiHu comes with a plot script that can plot 1D and 2D results (also over time) directly from these files (simply execute `plot` in the shell, the script is located under ``opendihu/scripts`` which is on ``PATH``).
  * **Python callback:** 
  
    You can provide a callback function in the input config file that will be called regularly with the previously described python dict.

  * **VTK/Paraview:** 
  
   Files with extensions `*.vtr`, `*.vts` and `*.vtu` are generated for the three mesh types, respectively and can be visualized with `ParaView <https://www.paraview.org/>`_. This output method is preferred for 3D results. ParaView provides extensive tools for manipulating the visualisation and post)-processing the results.

  * **Exfiles:**
  
    Every mesh type (not just `UnstructuredDeformableOfDimension<D>`) can be exported to pairs of `exelem` and `exnode` files. A `com` script is created that loads all generated exfiles (with correct offsets) and can be directly visualized using `cmgui <https://physiomeproject.org/software/opencmiss/cmgui/download>`_. 

  * **ADIOS native files:**
    
    The `ADIOS library <https://csmd.ornl.gov/adios>`_ is used to write this type of files. To perform in-situ visualization with `MegaMol <https://megamol.org/>`_, the files can be directly written to RAM. This type of output is only suggested for experts from `VISUS <https://www.visus.uni-stuttgart.de/en>`_.

* CSV based log files that contain parameters, timing and numerical information will be written at the end of each simulation run using parallel file I/O. They can be used for runtime statistics.

Parallelism
----------------

* Distributed memory parallelism using MPI is implemented for structured meshes (`StructuredRegularFixedOfDimension<D>` and `StructuredDeformableOfDimension<D>`). 
* Input can be specified globally and locally distributed. 
  
  * The simple case is a global specification of the mesh. In parallel execution, every process picks only the information for its own subdomain.
  * The more advanced and parallel efficient way is to specify mesh information locally for each process. This approach is required for scenarios with large meshe and highly parallel execution (HPC), where providing the whole information to a single process is not feasible due to memory limits.
  
    In such case, a single settings file can be specified for every process. It contains if branches depending on the own process rank, to determine which data to provide for which process.
* Output of python files and VTK files is uses parallelism. Parallel output in this case results in seperate files for each process. The python plotting utility of OpenDiHu transparently handles parallel output files in the same way as serially created files. Paraview is able to work with the parallel VTK files.
* MPI I/O is used to write combined VTK output files, i.e. a single file per time step. This is needed on supercomputers when running with a high number of cores.
* The monodomain example has been successfully executed on 27,000 cores on Hazel Hen to simulate a biceps with a typical number of 270,000 fibers.
* Instruction level parallelism is enabled by suitable data structures. This holds especially in the CellML functionality, where multiple instances of the model are combined to enable Single-instruction-multiple-data type parallelism.
  