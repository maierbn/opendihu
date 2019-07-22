.. _installation:

Installation
=================
Opendihu uses existing open-source projects, like PETSc, Python, Easylogging++, etc. The installation of opendihu has to provide all these packages, too. 
A `scons <https://scons.org/>`_ based build system is included that automatically downloads and builds all needed dependencies. 
It was successfully tested on Ubuntu 16.04, 18.04 and Debian as well as on the supercomputer Hazel Hen. It should work on other Linux distributions as well. If something fails, usually minor adjustments in the configuration solve the problem.

For users that only want to quickly check the functionality without a lengthy installation process, we provide a docker image of opendihu. This serves all the available functionality, except that parallel execution in docker containers is generally hardly possible. Because this is key to efficiently computing simulations, we recommend the native installation.

Using docker
----------------
Using the docker image you can use the framework directly without having to build and install any dependencies. The disadvantage is that you only have a shell and can't plot anything.

First install `docker <https://docs.docker.com/install/linux/docker-ce/ubuntu/>`_, following the instructions on the website. Then, start a shell inside the container with the following command:

.. code-block:: bash

  docker run -it maierbn/opendihu:latest bash

Inside the container, run `cd opendihu`, `git pull` and `make` to get and compile the latest code.

Native installation
----------------------
In order to use the code for development or for more efficient runs, it is necessary to clone the repository locally and build and install the framework including all dependencies.

.. code-block:: bash

  git clone https://github.com/maierbn/opendihu

There are several branches. The `develop` branch contains a recent version and is more-or-less stable. The `stable` branch is always stable but does not contain the latest developments. There are also multiple feature branches.

There is one `release <https://github.com/maierbn/opendihu/releases>`_ so far: version 1.0 from 15.04.2019. 

Prerequisites
^^^^^^^^^^^^^^

As a prerequisite, on a blank machine with ubuntu (tested on 16.04 and 18.04) you need to install the following packages.

.. code-block:: bash

  # Install prerequisites
  sudo apt-get update && \
  sudo apt-get install -y libopenmpi-dev libx11-* python2.7 git apt-utils make software-properties-common zlib1g-dev cmake libssl-dev bison flex

Because we use C++14, GCC version 5 or higher is required including the gfortran compiler. Ubuntu 16.04 has GCC 4 as default compiler chain, so you need to update to GCC 5 as follows. For Ubuntu 18.04 and later, this step is not necessary.

.. code-block:: bash

  # Install GCC5 toolchain
  sudo add-apt-repository ppa:ubuntu-toolchain-r/test && \
  sudo apt-get update && \
  sudo apt-get install -y gcc-5 g++-5 gfortran-5 && \
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-5

Make sure that the `gfortran` compiler is installed as well:

.. code-block:: bash

  sudo apt-get install gfortran

The `scons` build system needs python2.7. Make sure that the command `python2.7` starts a python 2.7 shell. If not, you probably have to create the following symbolic link:

.. code-block:: bash

  # link python executable
  ln -s python2.7 /usr/bin/python

All other needed dependencies will be handled by the `scons` build system. For each dependency you can either specify the path of its installation, if the dependency package is already installed on your system. Or you can not specify anything and let the build system download, build and install the dependencies on its own.
Note that `python3` with `numpy`, ``scipy`` and `matplotlib` is such a dependency. Opendihu will download and install `python3` including these packages.

Build 
^^^^^^^^^^^

The installation procedure can be started by the command `scons BUILD_TYPE=debug` for debug build or simply `scons` for release build. For convenience, there is also a `Makefile` that builds debug followed by release mode. The recommended way for the first installation is to simply execute

.. code-block:: bash

  make

There is also `make debug` and `make release` that just call `scons` with the respective build type and thus building the framework `debug` or `release` mode.

Instead of using the `Makefile` you can call `scons` yourself.

.. _installation_aliases:

Building with scons
^^^^^^^^^^^^^^^^^^^^^^^^

In order to build examples you need to use `scons`. The opendihu library can either be build using `scons` or using the `Makefile`, which again simply calls scons.

So you can either install scons on your system or use the `scons` program, that is packaged with opendihu. This is located under `dependencies/scons/scons.py`. It needs to be run with python 2.7 (not python3). 

It is advisable to define a bash alias for this scons command.
If you like, you can copy the following aliases to your `~/.bashrc` or `~/.bash_aliases` file:

.. code-block:: bash

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

Then simply execute ``sd`` to build in debug or ``s`` to build in release mode. Other options are ``sdd`` to build an example in debug mode from within the `build_debug` directory or analogously ``srr`` for release mode.

If you have called `make` and the framework compiled after some hours (green text), you were successful. Go on and build some examples (See next page, :doc:`getting_started`).
If not, read on, to find out what you need to configure in your case.

Configuring the build
^^^^^^^^^^^^^^^^^^^^^^^^

Configuration settings have to be provided in the python script `user-variables.scons.py`.
For every dependency package there are variables like

.. code-block:: bash

  #PETSC_DOWNLOAD=True
  #PETSC_DIR="~/petsc/install"

(Note, `#` means commented out here, because you shouldn't specify both lines at once). The first line would instruct the build system to download and build the package, in this case PETSc. The second line would provide the path to an already existing installation on the system, which would then be used. Thus, specify either of those. 

There are similar options for all packages. You can read about more possibilities in the header of the `user-variables.scons.py` file. 

There are required dependencies, which need to be present in order for opendihu to work, and optional dependencies:

============================================================  ========  ===================================================================================
 Package                                                      Required    Description
============================================================  ========  ===================================================================================
`MPI`                                                             yes     | *Message Passing Interface*, used for data transfer between
                                                                          | processes. This should be your system MPI. If you let 
                                                                          | opendihu install it for you, `OpenMPI <https://www.open-mpi.org/>`_ 
                                                                          | will be chosen.
`LAPACK`, `BLAS`                                                  yes     | Parallel linear algebra functions, this is a prerequisite 
                                                                          | to *PETSc*. Opendihu will install `OpenBLAS <https://github.com/xianyi/OpenBLAS/wiki>`_
`PETSc <https://www.mcs.anl.gov/petsc/>`_                         yes     | Low-level data structures and solvers, see their `website <https://www.mcs.anl.gov/petsc/>`_
                                                                          | for more details.
`Python3`                                                         yes     | The `python interpreter <https://www.python.org/>`_, 
                                                                          | version 3.6.5. We need the development header and source 
                                                                          | files, therefore it is recommended to let opendihu build 
                                                                          | python for you, even if your system has python installed.
`pythonPackages`                                                  yes     | This is a custom collection of python packages for the
                                                                          | python 3 interpreter, which is later available in the
                                                                          | python configuration scripts. It consists of 
                                                                          | `numpy matplotlib scipy numpy-stl svg.path triangle`.
`Easylogging++ <https://github.com/zuhd-org/easyloggingpp>`_      yes     | The used logging library. By default, logs are created 
                                                                          | in `/tmp/logs/` and output to the standard output.
`Base64 <https://github.com/tkislan/base64>`_                     yes     | An encoding standard and library that is used to create
                                                                          | binary VTK output files that can be viewed in Paraview.
                                                                          | Base64 encoded data is ASCII characters, the size is 4/3
                                                                          | of the raw binary data. The advantage is that despite 
                                                                          | being packed, it can be embedded in human-readable `XML`
                                                                          | files, which is the concept of VTK files.
`googletest <https://github.com/google/googletest>`_              no      | A testing framework, used for unit tests. Opendihu
                                                                          | compiles also without unit tests, but it is recommended 
                                                                          | to have them, especially when developing within the core.
`SEMT <https://github.com/maierbn/semt>`_                         no      | This is a small C++ symbolic differentiation toolbox 
                                                                          | that will be used for nonlinear solid mechanics, to 
                                                                          | derive material laws.
`ADIOS2 <https://adios2.readthedocs.io/en/latest>`_               no      | Binary output file format and library, parallely 
                                                                          | efficient and self-descriptive. This only installs, 
                                                                          | if you have a very recent version of `cmake`. It is no
                                                                          | problem, if this fails to install as most users won't 
                                                                          | need it. It is needed for interfacing `MegaMol`.
`MegaMol <https://megamol.org/>`_                                 no      | The parallel visualization framework developed at VISUS,
                                                                          | Uni Stuttgart. This installs the official version. To 
                                                                          | interface with opendihu, you would need a version that 
                                                                          | is not yet released. Therefore it is fine, if this is
                                                                          | not installed.
============================================================  ========  =================================================================================== 

It is recommended to not let the build system download and build `MPI`, instead you should use your local MPI installation. 

On Ubuntu systems, the system MPI directory should already be set correctly by the default value in `user-variables.scons.py`. Now run `make` to see, if MPI will be found.

If the MPI location is not detected automatically, you have to specifiy the path. Find out in which path on your system MPI is installed. The required directory contains a `lib` and an `include` subdirectory. It may be located at `/usr/lib/openmpi`, `/usr/lib/mpich`, `/usr/lib/x86_64-linux-gnu/openmpi` or similar. Set this path in `user-variables.scons.py` as value of the variable `MPI_DIR`.

When running `make`, `make debug` or `make release`, the dependencies will be downloaded and installed, and consequently, debug or release target will be build. The installation of dependencies can take several hours. The compilation afterwards completes in several minutes.

Troubleshooting
^^^^^^^^^^^^^^^^^^

If something fails during the installation, read the `config.log` file, which will be created. It contains information about the build process.

The dependencies that were installed successfully will be detected the next time and not installed again. You can force to rebuild selected packages by the `..._REBUILD` option, e.g.

.. code-block:: bash

  scons PETSC_REBUILD=True

to rebuild petsc, even if it was already detected. The same options that can be specified in the `user-variables.scons.py` file can also be given like this on the command line.

To also download the package and then install it again, use the `..._REDOWNLOAD` option, like

.. code-block:: bash

  scons PETSC_REDOWNLOAD=True

Sometimes it also helps to delete the folder of a package in the `dependencies` subdirectory and retry the installation. 

If a dependency fails to install, you can try to install it manually on your own. The commands that are used by the `scons` build system are logged in the `config.log` file.

If you want to change the build system to update the commands that are executed for installing a specific dependency, have a look at the directory `opendihu/dependencies/scons-config/sconsconfig/packages`. It contains the source for the build system. The main implementation is in `Package.py`, all other classes inherit from this class. Usually you find the file that is named like the dependency, e.g., `LAPACK.py` for Lapack or `PETSc.py` for PETSc.

If you change something here, you need to rebuild the python `egg` file of `scons-config`:

.. code-block:: bash

  cd <your-opendihu-path>
  cd dependencies/scons-config
  . install_manually.sh

Then, rerun the installation from the `opendihu` directory with `scons`.

