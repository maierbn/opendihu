[![Build Status](https://travis-ci.com/maierbn/opendihu.svg?branch=develop)](https://travis-ci.com/maierbn/opendihu)
[![CodeFactor](https://www.codefactor.io/repository/github/maierbn/opendihu/badge/develop)](https://www.codefactor.io/repository/github/maierbn/opendihu/overview/develop)

# Overview
The working title of this software framework is "opendihu" - from the project name "Digital Human". It serves as a code base to solve static and dynamic problems, where the Finite Element Method is used for spatial discretization. Due to its modular nature it is supposed to be adaptible for future problems.

This branch is for evaluation of different quadrature schemes for solid mechanics problems.

# Installation
Linux is required with python2.7, gcc version 5 or higher, gfortran (only needed when LAPACK is build)

On a blank machine with ubuntu you should install the following. Note that GCC 5 or higher is required.
```sudo apt-get install git libopenmpi-dev libx11-* python2.7`
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc-5 g++-5 gfortran-5
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-5```
Depending on your system you might already have all of that.

* Find out in which path on your system MPI is installed. A directory that contains a `lib` and an `include` subdirectory is needed. Often is is named like `/usr/lib/openmpi` or `/usr/lib/mpich`. 
* Set this path in `user-variables.scons.py` at line 20. You can also set paths for other packages, if they are already installed. This would reduce build time, it is, however, not required.
* Type `make debug` to build the debug version or `make release` to build the release version. This will at first download and install several needed dependencies. It takes a while, sometimes some components fail to install.
* You find some examples in the `examples` directory. Also the system tests under `testing/system_testing/tests` might be useful.
