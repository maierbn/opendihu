# Overview
The working title of this software framework is "opendihu" - from the project name "Digital Human". It serves as a code base to solve static and dynamic problems, where the Finite Element Method is used for spatial discretization. Due to its modular nature it is supposed to be adaptible for future problems.

This branch is for evaluation of different quadrature schemes for solid mechanics problems.

# Installation
Linux is required with python2.7.

* Find out in which path on your system MPI is installed. A directory that contains a `lib` and an `include` subdirectory is needed. Often is is named like `/usr/lib/openmpi` or `/usr/lib/mpich`. 
* Set this path in `user-variables.scons.py` at line 20
* Type `make debug` to build the debug version or `make release` to build the release version. This will at first download and install several needed dependencies. It takes up to 15 min, sometimes some components fail to install.
* You find some examples in the `examples` directory. The quadrature example is in the subdirectory `quadrature`. In that directory type `make` to build the example.
