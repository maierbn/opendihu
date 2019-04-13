[![Build Status](https://travis-ci.com/maierbn/opendihu.svg?branch=develop)](https://travis-ci.com/maierbn/opendihu)
[![CodeFactor](https://www.codefactor.io/repository/github/maierbn/opendihu/badge/develop)](https://www.codefactor.io/repository/github/maierbn/opendihu/overview/develop)

# Overview
Opendihu is a software framework that solves static and dynamic multi-physics problems, spatially discretized in 1D, 2D and 3D by the finite element method.
Our core design goals are usability, performance and extensibility.

Good usability refers to the process of configuring a given simulation with the python interface, running variants of the simulation and producing output in helpful file formats. Because of the contained python interpreter reconfiguration is possible at runtime. 

The performance goal is satisfied within the C++ core implementation. The data structures avoid expensive data copy and allow for vectorization. All structured grid functionality is designed for parallel execution. We closely build on the parallely efficient [PETSc](https://www.mcs.anl.gov/petsc/) library without overhead. The code was successfully run on 27,000 cores on the supercomputer "Hazel Hen" in Stuttgart.

The framework is extensible to future models. It provides infrastructure to store and manipulate scalar and vector fields, handle input and output, assemble system matrices and interface various PETSc solvers. For example, there is support for [CellML](https://www.cellml.org/), a format for interchanging systems of ordinary equations. 
The framework and its applications are constantly extended. However, there are stable releases.

# Installation
Opendihu is built upon existing, successful open-source projects, like PETSc, Python, Easylogging++, etc. This means that an installation has to provide all these packages, too. 
We include a build system that automatically downloads and builds all needed dependencies. It was successfully tested on Ubuntu 16.04, 18.04 and Debian as well as on the supercomputer Hazel Hen. It should work on other Linux distributions as well. If something fails, usually minor adjustments in the configuration solve the problem.

For users that only want to quickly check the functionality without a lengthy installation process, we provide a docker image of opendihu. This serves all the available functionality, except that parallel execution in docker containers is generally hardly possible. Because this is key to efficiently computing simulations, we recommend the native installation.

## Using docker
Using the docker image you can use the framework directly without having to build and install any dependencies. First install [docker](https://docs.docker.com/install/linux/docker-ce/ubuntu/). Then, start a shell inside the container by the following command:
```
docker run -it maierbn/opendihu_system_testing:latest bash
```
In the container, run `git pull` and `make` to get and compile the latest code.

## Native installation
In order to use the code for development or for more efficient runs, it is necessary to clone the present repository locally and build and install the framework including all dependencies.
(If you usually don't read instructions, clone and run `make` in the top level directory and see what happens. If there are error messages, look into the log file `config.log`.)

As a prerequisite, on a blank machine with ubuntu (tested on 16.04 and 18.04) you need to install the following packages.

```
  # Install prerequisites
  sudo apt-get update && \
  sudo apt-get install -y libopenmpi-dev libx11-* python2.7 git apt-utils make software-properties-common zlib1g-dev cmake libssl-dev bison flex
```

Because we use C++14, GCC version 5 or higher is required including the gfortran compiler. Ubuntu 16.04 has GCC 4 as default compiler chain, so you need to update to GCC 5 as follows. For Ubuntu 18.04 and later, this step is not necessary.

```
  # Install GCC5 toolchain
  sudo add-apt-repository ppa:ubuntu-toolchain-r/test && \
  sudo apt-get update && \
  sudo apt-get install -y gcc-5 g++-5 gfortran-5 && \
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-5
```

Make sure that the `gfortran` compiler is also installed:

```
  sudo apt-get install gfortran
```
The scons build system needs python2.7. Make sure that the command `python2.7` starts a python 2.7 shell. If not, you probably have to create the following symbolic link:

```
  # link python executable
  ln -s python2.7 /usr/bin/python
```

All other needed dependencies are handled by the `scons` build system. For each dependency you can either specify the path of its installation, if it is already present on your system or do not specify anything and let the build system download, build and install it on its own.
Note that python3 with numpy, scipy and matplotlib is such a dependency. Opendihu will download and install python3 and these packages.

The installation procedure is later started by the command `scons BUILD_TYPE=debug` for debug build or simply `scons` for release build. There is also a `Makefile` with the default target building debug mode and then release mode. Thus, the recommended way of the first installation is by running
```
make
```
There is also `make debug` and `make release` that just call `scons` with the respective build type.
But before doing this, you should read the following instructions about configuration.

Configuration settings have to be provided in the python script `user-variables.scons.py`.
For every dependency package there are variables like
```
#PETSC_DOWNLOAD=True
#PETSC_DIR="~/petsc/install"
```
(Note, `#` means commented out here, because you shouldn't specify both lines at once). The first line would instruct the build system to download and build the package, in this case PETSc. The second line would provide the path to an already existing installation on the system, which would then be used. Thus, specify either of those. 

There are similar options for all packages. You can read about more possibilities in the header of the `user-variables.scons.py` file. 

The required dependencies that need to be present in order for opendihu to work, are:

| Package | Required | Description |
| --------|----------|------------ |
| `MPI`   | yes      | *Message Passing Interface*, used for data transfer between processes. This should be your system MPI, if you let opendihu install it for you, [*OpenMPI*](https://www.open-mpi.org/) will chosen. |
| `LAPACK`,`BLAS` | yes | Parallel linear algebra functions, this is a prerequisite to *PETSc*. Opendihu will install [*OpenBLAS*](https://github.com/xianyi/OpenBLAS/wiki). |
| `PETSc` | yes      | Low-level data structures and solvers, see their [website](https://www.mcs.anl.gov/petsc/) for more details. |
| `Python3` | yes    | The [python interpreter](https://www.python.org/), version 3.6.5. We need the development header and source files, therefore it is recommended to let opendihu build python for you, even if your system has python installed. |
| `pythonPackages` | yes  | This is a custom collection of python packages for the python 3 interpreter, that are later available in the python configuration scripts. It consists of `numpy matplotlib scipy numpy-stl svg.path triangle`. |
| [`Easylogging++`](https://github.com/zuhd-org/easyloggingpp) | yes | The used logging library. By default, logs are created in `/tmp/logs/` and output to the standard output. |
| [`Base64`](https://github.com/tkislan/base64) | yes   | An encoding standard and library that is used to create binary VTK output files that can be viewed in Paraview. Base64 encoded data is ASCII characters, the size is 4/3 of the raw binary data. The advantage is that is packed and can be embedded in human-readable `XML` files, which is the concept of VTK files. |
| [`googletest`](https://github.com/google/googletest) | no   | A testing framework, used for unit tests. Opendihu obviously compiles also without unit tests, but it is recomme*nded to have them, especially when developing within the core. |
| [`SEMT`](https://github.com/maierbn/semt) | no     | This is a small C++ symbolic differentiation toolbox that will be used for nonlinear solid mechanics, to derive material laws. |
| [`ADIOS2`](https://adios2.readthedocs.io/en/latest) | no | Binary output file format and library, parallely efficient and self-descriptive. This only installs, if you have a very recent version of `cmake`. It is no problem, if this fails to install as most users won't need it. It is needed for interfacing `MegaMol`. |
| [`MegaMol`](https://megamol.org/) | no    | The parallel visualization framework developed at VISUS, Uni Stuttgart. This installs the official version. To interface with opendihu, you would need a version that is not yet released. Therefore it is fine, if this is not installed. |

It is recommended to not let the build system download and build `MPI`, instead you should use your local MPI installation. 
* On Ubuntu systems, the system MPI directory should already be set correctly by the default value in `user-variables.scons.py`. Now run `make` to see, if MPI will be found.
* If the MPI location is not detected automatically, you have to specifiy the path. Find out in which path on your system MPI is installed. The required directory contains a `lib` and an `include` subdirectory. It may be located at `/usr/lib/openmpi`, `/usr/lib/mpich`, `/usr/lib/x86_64-linux-gnu/openmpi` or similar.
* Set this path in `user-variables.scons.py` as value of the variable `MPI_DIR`.

When running `make`, `make debug` or `make release`, the dependencies will be downloaded and installed, and then the debug or release target will be build.
The installation of dependencies takes several hours. The compilation afterwards completes in some minutes. 
* If something fails, read the `config.log` file, which will be created. It contains information about the build process. Sometimes it helps to delete the folder of a package in the `dependencies` subdirectory and retry the installation. 
* The dependencies that were installed successfully will be detected the next time and not installed again. You can force to rebuild selected packages by the `..._REBUILD` option, e.g.
  ```
  scons PETSC_REBUILD=True
  ```
to rebuild petsc, even if it was already detected. The same options that can be specified in the `user-variables.scons.py` file can also be given like this on the command line.
* To also download the package and then install it again, use the `..._REDOWNLOAD` option, like
  ```
  scons PETSC_REDOWNLOAD=True
  ```
* If you call scons directly (instead of using the `make` wrapper), you can either install it on your system or use the scons, that comes with opendihu. It is located in `dependencies/scons/scons.py`. It needs to be run with python 2.7 (not python3). Then it might be useful to define an alias. If you like, you can copy the following to your `~/.bashrc` or `~/.bash_aliases` file:
```
alias scons='<your path>/opendihu/dependencies/scons/scons.py'
alias s='scons'
alias sd='scons BUILD_TYPE=d'
alias sdd='cd .. && scons BUILD_TYPE=d; cd -'
alias sddn='cd .. && scons BUILD_TYPE=d no_tests=yes no_examples=yes; cd -'
alias sdn='scons BUILD_TYPE=d no_tests=yes no_examples=yes'
alias srn='scons BUILD_TYPE=r no_tests=yes no_examples=yes'
alias sr='scons BUILD_TYPE=r'
alias srr='cd .. && scons BUILD_TYPE=r; cd -'
alias sdr='scons BUILD_TYPE=rd'
alias srd='scons BUILD_TYPE=rd'
alias srdd='cd .. && scons BUILD_TYPE=rd; cd -'
```

## Getting started
* To get started you find some examples in the `examples` directory. Also the system tests under `testing/system_testing/tests` might be useful.
* To build an example, `cd` into a subdirectory under `examples`, e.g. `examples/laplace/laplace2d`. In this directory run `scons`. 
  For this to work, you either need to install `scons` on your system (e.g. `sudo apt install scons` on ubuntu). Or you use the given `scons` in the `dependencies` directory (see above): 
```
  python2.7 ../../dependencies/scons/scons.py 
```
* To build the release target, use `scons` or `scons BUILD_TYPE=release` or `scons BUILD_TYPE=r`, to build the debug target, use `scons BUILD_TYPE=debug` or `scons BUILD_TYPE=d`.
* There will be executables created in the `build_debug` or `build_release` subdirectories. Change into one of these directories and run the program with a settings file as only argument: `./laplace_regular settings_lagrange_quadratic.py`.
* Output files in this example (and likewise in the other examples) will be created under the `out` subdirectory. If you look into `out` you'll find two files: `laplace.py` and `laplace.vtr`.
 
  The `*.vtr` files can be visualized using Paraview. The `*.py` files can be visualized using the plot script in `opendihu/scripts`. 
  It is useful to add the `scripts` directory to the `PATH` environment variable, e.g. by adding the line `export PATH=$PATH:<your-path>/opendihu/scripts` to your `.bashrc` file.
  Then you can run `plot.py <*.py-output files>` anywhere to visualize the output. Usually the shortcut `plot` instead of `plot.py` should also work. (Unless you have installed something different with the name `plot`).
  
  An often used command is therefore `plot out/*`. If arguments are ommited, i.e. `plot`, this is the same as `plot *.py`.
  In our example the command could be `plot out/laplace.py`.
* The source files are located in the `src` subdirectory. The files to be compiled are specified in `SConscript`. 
  There, you can, for example, comment out the examples that you don't want to compile everytime.
* Now, change into `src` (i.e. into `.../examples/laplace/laplace2d/src`) and open `laplace_regular.cpp`. This is the main file of the 2D Laplace model. As can be seen, the equation `Δu = 0` is discretized by the Finite Element Method on a structured regular grid of dimension 2, basis functions are Lagrange of order 2, and the Quadrature scheme is Gauss quadrature with 3 gauss points per dimension. 
  Now change to linear Lagrange basis functions, by changing `LagrangeOfOrder<2>` to `LagrangeOfOrder<1>`.
  Change into the `build_debug` directory (`cd ../build_debug`). 
  If you have set the aliases of Sec. 2, you can recompile with `sdd`. Otherwise go up one directory and run `scons BUILD_TYPE=d`. 
  Now, from the `build_debug` directory, run the new executable with 
```
./laplace_regular ../settings_lagrange_linear.py
```
  Plot the result with `plot out/*`.
* Test the parallel execution and run the same program with the same settings file on two processes:
```
mpirun -n 2 ./laplace_regular ../settings_lagrange_linear.py
```
  If you now look into the out directory (`ls -l out`), you'll see that two new files `laplace.0.py` and `laplace.1.py` were created from the two processes. The file `laplace.py` is still the old one from the single process.
  Now plot the new files, either `plot out/laplace.0.py out/laplace.1.py` or shorter `plot out/laplace.*.py`. The result looks the same.
  Check that the results from the serial and parallel are actually the same using the following helper script:
```
validate_parallel.py out/*
```
* The created python output files are human-readable (because `"binary":False` is set in the settings file). You can open them in an editor and see what they contain. There is also a script for formatted printing on the console:
```
catpy out/laplace.0.py
```
* With the current settings, also the Paraview files are human-readable. You can also open e.g. `out/laplace.vtr` in an editor. Also try loading the `.pvtr` file in Paraview. 
  For big files it is better to produce binary files.
  In the settings file `settings_lagrange_linear.py` change `"binary":False` to `"binary":True` in the output writers. Now if you run the program again you'll get binary files that can't be read in a text editor. However, the `plot`, `validate_parallel` and `catpy` utilities still work. 
* If you know `cmgui`, the visualization tool of OpenCMISS Zinc, you can also generate `exnode` and `exelem` output files for cmgui. Add the line
```
  {"format": "Exfile", "filename": "out/laplace"},
```
to the `"OutputWriter"` list in file `settings_lagrange_linear.py` (line 31). After running the program again, you get the output files `laplace.exelem`, `laplace.exnode` and `laplace.com` in the out directory. The `.com` file is a convienient perl script that sets up the visualization in cmgui (OpenCMISS Iron won't generate this for you.). Change into the out directory and simply run `cmgui laplace.com`. In the Scene Editor click on `/` and then the `surface` item. Under `data`, select `solution` as the field variable that will be shown in color. Now you can tilt the view in the Graphics window to see the solution.
* Now you know the basics how to run a simulation program. You can try to change parameters in the settings file, like number of elements (variables `m` and `n`), the `physicalExtent` or try to understand, how the Dirichlet boundary conditions were specified. Note, that because this example uses a `Mesh::StructuredRegularFixedOfDimension<2>` mesh (in the `cpp` source file), we can only have elements with quadratic shape, i.e. `physicalExtent` and `nElements` have to match. You can look into the `laplace_structured.cpp` example file, which uses a structured mesh, that can have different mesh width in `x` and `y` direction or even arbitrary node positions.
* The settings files use python syntax and are actually python scripts. This means you can execute any python code there, for example load your own custom geometry or input data files and set the options appropriately. The general documentation of the options is only given through the examples, so if you need to know how to specify certain options, look for an example files, that does it, or ask me.

# Documentation
Theory documentation can be found in the `doc/derivations/doc.pdf` document. Some developer hints for the core code can be found in `doc/documentation.rst`. Generally the C++ code contains a lot of comments and often it is useful to look directly into the code, to find out how something works. The real documentation is still in the developers, so ask if you need to know something.
The following functionality is currently implemented:

## Equations
Supported equation types are currently
* Laplace and Poisson
* Diffusion, istropoic and anisotropic
* Monodomain: Reaction-diffusion, where the reaction term is given by a CellML description
* Bidomain equation
* Multidomain equation
* Multi-physics combinations of the above
* Functionality to estimate and trace fibers from streamlines of a laplacian flow through a muscle volume

## Boundary and initial conditions
* Dirichlet and Neumann-type boundary conditions are supported

## Discretization and Mesh types
* For the time domain, there are explicit and implicit Euler and Heun method implemented, Crank-Nicolson, as well as a Godunov and Strang splitting scheme. 
* There are three implemented mesh types which each have variants for 1D, 2D and 3D.
  * `StructuredRegularFixedOfDimension<D>`: This is a rectilinear grid with fixed mesh width in all dimensions. To store a mesh, only the number of elements in the `D` dimensions and the mesh width or extent is needed. This is computationally fast and memory efficient. The stiffness and mass matrices for the laplace operator are assembled by precomputed stencils. This mesh type is mainly used for validation of the other mesh types and for examples that have an analytical solution.
  * `StructuredDeformableOfDimension<D>`: This structured mesh can have arbitrary node positions, and there is a fixed number of elements in each dimension. Internally there is no need to store extra adjacency information, only the node positions and number of elements in the dimensions is stored. This is also quite fast. This mesh can deform over time, therefore it is suitable for structured mechanics computations.
  * `UnstructuredDeformableOfDimension<D>`: This is the most general mesh type. The elements can arbitrarily share nodes. Supported are also *versions* like in OpenCMISS where there can be multiple degrees of freedom at a node, to model incontinuities. These meshes can be natively constructed from exfiles.
* Available basis functions are linear and quadratic Lagrange and cubic Hermite. Unit tests cover all combinations of mesh types in 1D,2D,3D with the basis function types, in serial and parallel execution.

# Input/Output
* The framework has a Python3 interpreter included and defines its own config format by a python dictionary. Input files of meshes can be either given by python lists of node positions, element-node mappings, number of elements etc. (depending on the mesh type). By this is it possible to parse further data sources in the input config file which will be parsed by the Python interpreter at the beginning of the program.
* Also `exnode` and `exelem` files can be parsed to construct a mesh of type `UnstructuredDeformableOfDimension<D>`.
* There are 5 output formats supported:
  * Python file: a data representation is exported as a python dictionary to a file. This can be configured to be human readable (ascii), or binary (using the python `pickle` package). This is the most memory efficient and easy to handle output format. There exists also a plot script that can plot 1D and 2D results (also over time) directly from these files.
  * Python callback: You can provide a callback function in the input config file that will be called regularly with the previously described python dict.
  * Exfiles: Every mesh type (not just `UnstructuredDeformableOfDimension<D>`) can be exported to `exelem`, `exnode` files. Also a `com` script is created that loads all generated exfiles (with correct offsets) and can be directly visualized using `cmgui`. However some manual tweaking with the `com` file is sometimes required.
  * VTK/Paraview: Files with extensions `*.vtr`, `*.vts` and `*.vtu` can be generated for the three mesh types, respectively. These are the preferable output method for 3D results. Paraview provides extensive tools for manipulation of the visualisation and postprocessing of the results. Only a recent version of Paraview can directly show 1D meshes.
  * ADIOS native files, these can also be writen to RAM, to perform in-situ visualization with MegaMol
* CSV based Log files containing parameters, timing and numerical information are written at the end of each simulation run, using parallel file I/O.

# Parallelism
* Distributed memory parallelism using MPI is implemented for structured meshes (`StructuredRegularFixedOfDimension<D>` and `StructuredDeformableOfDimension<D>`). 
* The python output files as well as VTK output files are parallel, the Exfiles output is serial. The python plotting utility can handle the parallel output files.
* MPI I/O is used to write combined VTK output files, i.e. a single file per time step. This is needed on supercomputers when running with a high number of cores.
* The monodomain example has been successfully executed on 27,000 cores on Hazel Hen to simulate a biceps with a typical number of 270,000 fibers.
  
# Tests
* There are unit tests that run after each compilation (you can abort the compilation process after the library was created to skip the unit tests). 
  There are also system tests that run longer scenarios for various settings and for do some comparisons to analytical solutions. The list of system tests is to be extended, currenty it only includes Laplace and Diffusion examples (but for all combinations of ansatz functions and mesh types). The system tests compile latex slides and a pdf document containing test results and also images and movies of the test runs. It runs nighly on a local jenkins installation. 
