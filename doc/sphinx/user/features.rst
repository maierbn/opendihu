
Features
============

This page lists the functionality that is currently implemented. More detailed information is given on the :doc:`/settings` pages.

Equations
-----------
Supported equations are currently:

* Laplace and Poisson
* Diffusion, istropoic and anisotropic
* Monodomain: Reaction-diffusion, where the reaction term is given by a CellML description
* Bidomain equation
* Multidomain equation
* Multi-physics combinations of the above
* Functionality to estimate and trace fibers from streamlines of a laplacian flow through a muscle volume
* Static linear elasticity

Boundary and initial conditions
-------------------------------------

* Dirichlet and Neumann-type boundary conditions are supported.
* For transient problems, initial conditions have to be specified.

Discretization and Mesh types
------------------------------------
* For the time domain, there are explicit and implicit Euler and Heun method implemented, Crank-Nicolson, as well as a Godunov and Strang splitting scheme. 
  To learn more about the time stepping schemes, see :doc:`/settings/timestepping_schemes_ode`, the splitting schemes are explained at :doc:`/settings/splitting`.
* There are three implemented mesh types which each have variants for 1D, 2D and 3D.

  * `StructuredRegularFixedOfDimension<D>`: This is a rectilinear grid with fixed mesh width in all dimensions. To store a mesh, only the number of elements in the `D` dimensions and the mesh width or extent is needed. This is computationally fast and memory efficient. The stiffness and mass matrices for the laplace operator are assembled by precomputed stencils. This mesh type is mainly used for validation of the other mesh types and for examples that have an analytical solution.
  * `StructuredDeformableOfDimension<D>`: This structured mesh can have arbitrary node positions, and there is a fixed number of elements in each dimension. Internally there is no need to store extra adjacency information, only the node positions and number of elements in the dimensions is stored. This is also quite fast. This mesh can deform over time, therefore it is suitable for structured mechanics computations.
  * `UnstructuredDeformableOfDimension<D>`: This is the most general mesh type. The elements can arbitrarily share nodes. Supported are also *versions* like in OpenCMISS where there can be multiple degrees of freedom at a node, to model incontinuities. These meshes can be natively constructed from exfiles.
* Available basis functions are linear and quadratic Lagrange and cubic Hermite. Unit tests cover all combinations of mesh types in 1D,2D,3D with the basis function types, in serial and parallel execution.

More details on the meshes are given on the :doc:`/settings/mesh` page.

Input/Output
---------------

* The framework has a Python3 interpreter included and defines its own config format by a python dictionary. Input files of meshes can be either given by python lists of node positions, element-node mappings, number of elements etc. (depending on the mesh type). By this is it possible to parse further data sources in the input config file which will be parsed by the Python interpreter at the beginning of the program.
* Also `exnode` and `exelem` files can be parsed to construct a mesh of type `UnstructuredDeformableOfDimension<D>`.
* There are 5 output formats supported:

  * **Python file:** a data representation is exported as a python dictionary to a file. This can be configured to be human readable (ascii), or binary (using the python `pickle` package). This is the most memory efficient and easy to handle output format. There exists also a plot script that can plot 1D and 2D results (also over time) directly from these files.
  * **Python callback:** You can provide a callback function in the input config file that will be called regularly with the previously described python dict.
  * **Exfiles:** Every mesh type (not just `UnstructuredDeformableOfDimension<D>`) can be exported to `exelem`, `exnode` files. Also a `com` script is created that loads all generated exfiles (with correct offsets) and can be directly visualized using `cmgui`. 
  * **VTK/Paraview:** Files with extensions `*.vtr`, `*.vts` and `*.vtu` can be generated for the three mesh types, respectively. These are the preferable output method for 3D results. Paraview provides extensive tools for manipulation of the visualisation and postprocessing of the results. Only a recent version of Paraview can directly show 1D meshes.
  * **ADIOS native files:** these can also be writen to RAM, to perform in-situ visualization with MegaMol
* CSV based log files that contain parameters, timing and numerical information, will be written at the end of each simulation run using parallel file I/O.

Parallelism
----------------

* Distributed memory parallelism using MPI is implemented for structured meshes (`StructuredRegularFixedOfDimension<D>` and `StructuredDeformableOfDimension<D>`). 
* Input can either be specified globally, then every process picks the information it needs from the global settings file. It is also possible to only specify the local information for each process. This is needed for large scenarios at highly parallel execution, where providing the whole information to a single process is not feasible. The settings file can still be the same file for every process. It can contain switches on the own process number, which data to provide in the config.
* The python output files as well as VTK output files are parallel, the Exfiles output is serial. The python plotting utility can handle the parallel output files.
* MPI I/O is used to write combined VTK output files, i.e. a single file per time step. This is needed on supercomputers when running with a high number of cores.
* The monodomain example has been successfully executed on 27,000 cores on Hazel Hen to simulate a biceps with a typical number of 270,000 fibers.
* Instruction level parallelism is enabled by suitable data structures. This holds especially in the CellML functionality, where multiple instances of the model are combined to enable Single-instruction-multiple-data type paralelism.
  
Tests
---------

* There are unit tests that run after each compilation (you can abort the compilation process after the library was created to skip the unit tests).   
  The continuous integration service Travis CI automatically builds and executes all unit tests after each push to the repo. This takes around 30 min, if tests fail, the responsible developer is notified via e-mail.
* There are also system tests that run longer scenarios for various settings and for do some comparisons to analytical solutions. The list of system tests is to be extended, currenty it only includes Laplace and Diffusion examples (but for all combinations of ansatz functions and mesh types). The system tests compile latex slides and a pdf document containing test results and also images and movies of the test runs. It runs nighly on a local jenkins installation. 
