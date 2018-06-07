[![Build Status](https://travis-ci.com/maierbn/opendihu.svg?branch=develop)](https://travis-ci.com/maierbn/opendihu)
[![CodeFactor](https://www.codefactor.io/repository/github/maierbn/opendihu/badge/develop)](https://www.codefactor.io/repository/github/maierbn/opendihu/overview/develop)

# Overview
The working title of this software framework is "opendihu" - from the project name "Digital Human". It serves as a code base to solve static and dynamic problems, where the Finite Element Method is used for spatial discretization. Due to its modular nature it is supposed to be adaptible for future problems.

# Installation
Linux is required with python2.7, gcc version 5 or higher, gfortran (only needed when LAPACK is build)

On a blank machine with ubuntu you should install the following. Note that GCC 5 or higher is required.
```
  sudo apt-get install git libopenmpi-dev libx11-* python2.7
  sudo add-apt-repository ppa:ubuntu-toolchain-r/test
  sudo apt-get update
  sudo apt-get install gcc-5 g++-5 gfortran-5
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-5
```
Depending on your system you might already have all of that.

* Find out in which path on your system MPI is installed. A directory that contains a `lib` and an `include` subdirectory is needed. Often is is named like `/usr/lib/openmpi` or `/usr/lib/mpich`. 
* Set this path in `user-variables.scons.py` at line 20. You can also set paths for other packages, if they are already installed. This would reduce build time, it is, however, not required.
* Type `make debug` to build the debug version or `make release` to build the release version. This will at first download and install several needed dependencies. It takes a while, sometimes some components fail to install.
* You find some examples in the `examples` directory. Also the system tests under `testing/system_testing/tests` might be useful.

# Documentation
Currently some theory documentation can be found in the `doc/derivations/doc.pdf` document. 
You can find out about how the framework is instantiated by looking at the examples in the `example` directory and the system tests in the `testing/system_testing/tests` directory. Where the examples sometimes are out-of-date, the system_tests run regularly and should normally work. 

The following functionality is currently implemented:
## Equations
Supported equation types are currently
* Laplace and Poisson
* Diffusion, istropoic and anisotropic
* Reaction-diffusion by the CellML adapter
* The nonlinear solid/structural mechanics equations are not fully functional at the moment. They work for some test scenarios but not for others. Also there is no dynamic elasticity yet.

## Discretization and Mesh types
* The finite element is the only available spatial discretization method.
* For the time domain, there are explicit Euler and Heun implemented, as well as a Godunov and Strang splitting scheme. The focus is currently not on numerical investigations. It should not be too hard to implement further schemes, like Crank-Nicholson.
* Ansatz functions are linear and quadratic Lagrange and cubic Hermite. The scenarios are essentially testing with all three types. Further ansatz function types should be easily extendable (like cubic Lagrange or different Hermite), but they're probably not really needed.
* There are three implemented mesh types which each have variants for 1D, 2D and 3D.
  * `StructuredRegularFixedOfDimension<D>`: This is a rectilinear grid with fixed mesh width in all dimensions. To store a mesh only the number of elements in the `D` dimensions and the mesh width is needed. This is computationally fast and memory efficient. The stiffness matrix for the laplace operator is assembled by precomputed stencils. This mesh type is mainly used for validation of the other mesh types and for some easy examples that have an analytical solution.
  * `StructuredDeformableOfDimension<D>`: This structured mesh can have arbitrary node positions, but there is a fixed number of elements in each dimension. Internally there is no need to store extra adjacency information, only the node positions and number of elements in the dimensions is stored. This is also quite fast. This mesh can deform over time, therefore it is suitable for structured mechanics computations.
  * `UnstructuredDeformableOfDimension<D>`: This is the most general mesh type. The elements can arbitrarily share nodes. Supported are also *versions* like in OpenCMISS where there can be multiple degrees of freedom at a node, to model incontinuities. These meshes can be natively constructed from exfiles.

# Input/Output
* The framework has a Python3 interpreter included and defines its own config format by a python dictionary. Input files of meshes can be either given by python lists of node positions, element-node mappings, number of elements etc. (depending on the mesh type). By this is it possible to parse further data sources in the input config file which will be parsed by the Python interpreter at the beginning of the program.
* Also `exnode` and `exelem` files can be parsed to construct a mesh of type `UnstructuredDeformableOfDimension<D>`.
* There are 4 output formats supported:
  * Python file: a data representation is exported as a python dictionary to a file. This can be configured to be human readable (ascii), or binary (using the python `pickle` package). This is the most memory efficient and easy to handle output format. There exists also a plot script that can plot 1D and 2D results (also over time) directly from these files.
  * Python callback: You can provide a callback function in the input config file that will be called regularly with the previously described python dict.
  * Exfiles: Every mesh type (not just `UnstructuredDeformableOfDimension<D>`) can be exported to `exelem`, `exnode` files. Also a `com` script is created that loads all generated exfiles (with correct offsets) and can be directly visualized using `cmgui`. However some manual tweaking with the `com` file is sometimes required.
  * Paraview: Files with extensions `*.vtr`, `*.vts` and `*.vtu` can be generated for the three mesh types, respectively. These are the preferable output method for 3D results. Paraview provides extensive tools for manipulation of the visualisation and postprocessing of the results. Only a recent version of Paraview can directly show 1D meshes.
  
# Tests
* There are unit tests that run after each compile (you can abort the compilation process after the library was created to skip the unit tests). There are also system tests that run real scenarios for multiple settings and for some compare to analytical solutions. The list of system tests is to be extended, currenty it only includes Laplace and Diffusion examples (but for all combinations of ansatz functions and mesh types). The system tests compile latex slides and a pdf document containing test results and also images and movies of the test runs. It runs nighly on a local jenkins installation.  
  
# Bug reporting
* If you find bugs, you can write me an e-mail, or even better you have a look at it and try to fix them on your own (make an own branch). Currently there are new features being added or internal code structure changes quite frequently, which sometimes introduces new bugs and fixes existing ones. There is some safety from the unit tests and system tests.

# Work list
* A major development step will be the parallelization using MPI/PETSc. For this the data structures are already somehow prepared, but some effort needs to be put into the dof mappings. (shared memory parallelisation using OpenMP is included already at some points)
* Also the continuum mechanics problem has to be addressed, still. (Probably will be before the parallelisation).
* The current work is concerned with multiple fibres on a realistic biceps brachii geometry.
