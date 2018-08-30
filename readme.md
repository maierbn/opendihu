[![Build Status](https://travis-ci.com/maierbn/opendihu.svg?branch=develop)](https://travis-ci.com/maierbn/opendihu)
[![CodeFactor](https://www.codefactor.io/repository/github/maierbn/opendihu/badge/develop)](https://www.codefactor.io/repository/github/maierbn/opendihu/overview/develop)

# Overview
The working title of this software framework is "opendihu" - from the project name "Digital Human". 
It serves as a code base to solve static and dynamic problems, where the Finite Element Method is used for spatial discretization. 
Due to its modular nature it is supposed to be adaptible for future problems.

# Installation
If you're impatient, type `make` in the top level directory and see what happens. If there are error messages, look into the log file `config.log`. 

But please read the following instructions first.

On a blank machine with ubuntu (tested on 16.04 and 18.04) you need to install the following packages.

```
    # Install prerequisites
    sudo apt-get update && \
    sudo apt-get install -y libopenmpi-dev libx11-* python2.7 git apt-utils make software-properties-common zlib1g-dev cmake libssl-dev
```

GCC version 5 or higher is required including the gfortran compiler. Ubuntu 16.04 has GCC 4 as standard compiler chain, so you need to update to GCC 5 as follows:

```
  # Install GCC5 toolchain
  sudo add-apt-repository ppa:ubuntu-toolchain-r/test && \
  sudo apt-get update && \
  sudo apt-get install -y gcc-5 g++-5 gfortran-5 && \
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-5
```

On Ubuntu 18.04 there is GCC 7 available, just install gfortran if don't have it already:

```
  sudo apt-get install gfortran
```
The scons build system needs python2.7. Make sure that the command `python2.7` starts a python shell. If not, create the following symlink:

```
  # link python executable
  ln -s python2.7 /usr/bin/python
```

All other needed dependencies are handled by the `scons` build system. For each dependency you can either specify the location of its installation or do nothing and let the build system download, build and install it for you.
Note that python3 with numpy, scipy and matplotlib is such a dependency. It may download and install python3 on its own.

It is recommended to not let the build system download and build `MPI`, instead specify the location of your local MPI installation.
* Find out in which path on your system MPI is installed. We need the directory that contains a `lib` and an `include` subdirectory. 
  It may be named like `/usr/lib/openmpi` or `/usr/lib/mpich` (e.g. on Ubuntu 16.04) or `/usr/lib/x86_64-linux-gnu/openmpi` (on Ubuntu 18.04).
* Set this path in `user-variables.scons.py` at line ~20. You don't have to do anything for Ubuntu 16.04 or Ubuntu 18.04, because it should already be done.
* In this file you can also set paths for other packages, if they are already installed. This would reduce build time, it is, however, not required.
* Type `make debug` or `make release` to build the debug or release target. If you run `make` without targets, it will build the debug target followed by the release target.
  At first it will download and install several needed dependencies. It takes several hours. If some components fail to install you can inspect the log file `config.log`. If only an optional component failed, opendihu may still work.
  The required dependencies are:
```
  MPI         (OpenMPI will be installed by default)
  LAPACK/BLAS (OpenBLAS will be installed by default)
  PETSc
  Python 3
  Base64
  Googletest
  SEMT
  Easylogging++
```
There are optional dependencies that allow compilation of opendihu. But they may be needed for unit tests. The following are the optional dependencies:
```
  bzip2
  Cython
  NumpyC
  Scipy
  Matplotlib
```
* To get started you find some examples in the `examples` directory. Also the system tests under `testing/system_testing/tests` might be useful.

# Documentation
The theory documentation can be found in the `doc/derivations/doc.pdf` document. Documentation concerning the code can be found in `doc/documentation.rst`.
The best way to learn how the framework is used is by looking at the examples in the `example` directory and the system tests in the `testing/system_testing/tests` directory. 

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
  * VTK/Paraview: Files with extensions `*.vtr`, `*.vts` and `*.vtu` can be generated for the three mesh types, respectively. These are the preferable output method for 3D results. Paraview provides extensive tools for manipulation of the visualisation and postprocessing of the results. Only a recent version of Paraview can directly show 1D meshes.

# Parallelism
* Distributed memory parallelism using MPI is implemented for structured meshes (`StructuredRegularFixedOfDimension<D>` and `StructuredDeformableOfDimension<D>`). 
* The python output files as well as VTK output files are parallel, the Exfiles output is only serial. The python plotting utility can handle the parallel output files.
  
# Tests
* There are unit tests that run after each compile (you can abort the compilation process after the library was created to skip the unit tests). 
  There are also system tests that run real scenarios for multiple settings and for some compare to analytical solutions. The list of system tests is to be extended, currenty it only includes Laplace and Diffusion examples (but for all combinations of ansatz functions and mesh types). The system tests compile latex slides and a pdf document containing test results and also images and movies of the test runs. It runs nighly on a local jenkins installation. 
