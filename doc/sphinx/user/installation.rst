.. _installation:

Installation
=================
OpenDiHu uses existing open-source projects, like PETSc, Python, Easylogging++, etc. The installation of opendihu has to provide all these packages, too. 
A `scons <https://scons.org/>`_ based build system is included that automatically downloads and builds all needed dependencies. 
It was successfully tested on Ubuntu 16.04, 18.04 and 20.04 (also on the Windows subsystem for linux, WSL) and on Debian as well as on the supercomputers `Hazel Hen` and `Hawk`. 
It should work on other Linux distributions as well. If something fails, usually minor adjustments in the configuration solve the problem.

For users that only want to quickly check the functionality without a lengthy installation process, we provide a docker image of opendihu.
This serves all the available functionality, except that parallel execution in docker containers is generally hardly possible. 
Because this is key to efficiently computing simulations, we recommend the native installation.

Using docker
----------------
Using the docker image you can use the framework directly without having to build and install any dependencies. 
The disadvantage is that you only have a shell and can't plot anything. The docker image may also be slightly out-of-date.

First install `docker <https://docs.docker.com/install/linux/docker-ce/ubuntu/>`_, following the instructions on the website. Then, start a shell inside the container with the following command:

.. code-block:: bash

  docker run -it maierbn/opendihu:latest bash

Inside the container, run ``cd opendihu``, ``git pull`` and ``make`` to get and compile the latest code.

Native installation
----------------------
In order to use the code for development or for more efficient runs, it is necessary to clone the repository locally and build and install the framework including all dependencies.

.. code-block:: bash

  git clone https://github.com/maierbn/opendihu

There are several branches. The `develop` branch contains a recent version and is more-or-less stable. This is the recommended branch to start with.
The `stable` branch aims to be always stable but does not contain the latest developments. There are also multiple feature branches.
Which branches are being developed can be checked on the overview page of the continuous integration tool, `CircleCI <https://app.circleci.com/pipelines/github/maierbn/opendihu>`_.

From time to time there are `releases <https://github.com/maierbn/opendihu/releases>`_, e.g. there is version 1.0 from 15.04.2019 and newer ones.

Prerequisites
^^^^^^^^^^^^^^

On a blank computer with Ubuntu, the following packages should be installed.
On all servers of SGS, these packages are already available and no further steps are required.

.. code-block:: bash

  # Packages needed on Ubuntu 16.04 and above
  sudo apt-get update && \
  sudo apt-get install -y git make cmake apt-utils software-properties-common libopenmpi-dev libx11-* zlib1g-dev libssl-dev libffi-dev bison flex

  # Additional packages needed on Ubuntu 20.04
  sudo apt install python-is-python3      # This makes `python` equal to `python3` and is only needed for the installation of PETSc.
  sudo apt install python3-pip
  
  # (optional) If you want to use precice, install the following (requires ~300 MB):
  sudo apt install libboost-filesystem-dev libboost-log-dev libboost-program-options-dev libboost-system-dev libboost-thread-dev libboost-test-dev
  
  # (optional) To be able to build this documentation, install
  sudo pip3 install sphinx recommonmark sphinx_rtd_theme


Because we use C++14, the GCC compiler **version 7 or higher** is required, including the Fortran compiler, ``gfortran``.
The default GCC version on Ubuntu 16.04 was GCC 4 and updating was required. From Ubuntu 18.04 on, the default version is at least version 7.

Please check your current version by running:

.. code-block:: bash

  gcc --version
  gfortran --version

Only if the version is GCC 6 or lower, you need to update your GCC. This can be done as follows.

.. code-block:: bash

  # Install GCC9 toolchain
  sudo add-apt-repository ppa:ubuntu-toolchain-r/test && \
  sudo apt-get update && \
  sudo apt-get install -y gcc-9 g++-9 gfortran-9 && \
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-9 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-9

Make sure that the `gfortran` compiler is installed as well:

.. code-block:: bash

  sudo apt-get install gfortran

All other needed dependencies will be handled by the `scons` build system. 
For each dependency you can either specify the path of its installation if the dependency package is already installed on your system.
Or, if you don't do anything special, the build system downloads, builds and installs the dependencies on its own. This is the recommended way.

Note that one of these dependencies is a development version of `python3` with `numpy`, `scipy` and `matplotlib`. OpenDiHu will download and install `python3` including these packages regardless of an already existing python3 installation on your system.

Build 
^^^^^^^^^^^
The recommended way for the first installation is to change into the `opendihu` directory and simply execute

.. code-block:: bash

  make

Then, scons will download and install everything for a while. It runs the unit tests using 1, 2 and 6 processes. Then, it compiles all examples. As soon as the unit tests are being compiled, the installation has finished and you can abort the process. Or, you can wait for it to finish.

If some of the dependencies were not found this is not a problem, e.g. if precice fails, you'll not be able to use precice but everything else still works.

You can also execute `make release` to only build the release target. This is enough if you don't aim at developing the C++ code.

Three different targets are defined: `release`, `debug` and `releasewithdebuginfo`. In `release` target, the code will be optimized to run as fast as possible.
In `debug` target, compilation and execution will take more time. A lot of debugging information will be printed by the program to the console. This is the standard target to use during development.
The third target, `releasewithdebuginfo` enables optimizations, like the `release` target, but additionally includes the debugging output.

Analogous to ``make release``, there is also ``make debug`` to build the debug target and ``make release_without_tests`` or `make debug_without_tests` to exclude build unit tests (which are not required but take a lot of time).
To learn about more available make targets, read the `Makefile`.

Internally, ``make`` calls the build system, `scons`.
The installation procedure can also be started by the command `scons` for release build or `scons BUILD_TYPE=debug` for debug build. 
The ``make`` targets ``make release`` and ``make debug`` just call ``scons`` with the respective build type and thus building the framework `debug` or `release` mode.
Instead of using the `Makefile` you can also call ``scons`` yourself.

.. _installation_aliases:

Define environment variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order for some commands to work (e.g. the ``plot`` utility), you need to set the PATH variable to point to some directories of opendihu. 
This can be done by adding the following lines to your `~/.bashrc` script or `~/.bash_aliases` on Ubuntu.

.. code-block:: bash

  # set environment variables and PATH
  export OPENDIHU_HOME=~/opendihu         # replace this by the location for your installation
  export PATH=$PATH:$OPENDIHU_HOME/scripts
  export PATH=$PATH:$OPENDIHU_HOME/scripts/geometry_manipulation
  export PATH=$PATH:$OPENDIHU_HOME/scripts/file_manipulation

(Replace the `~/opendihu` with your own path).
Setting these variables is recommended but not required.

The `~/.bashrc` or `~/.bash_aliases` file will be executed whenever you start a new `bash` instance. 
In order for the variable assignments to take effect, either close and reopen the console window or source the file yourself, by executing ``. ~/.bashrc``.

Building with scons
^^^^^^^^^^^^^^^^^^^^^^^^

Opendihu consists of a `core` library that contains the main functionality and multiple examples, that each use the core library.
As mentioned, to build the opendihu core library either `make` can be used, or it is possible to use the build system `scons`.
In order to build examples there is no choice, you need to use `scons`.

To be able to use `scons`, you can either install the `scons` package on your system (``sudo apt install scons`` on Ubuntu)
or use the `scons` program, that is packaged with opendihu. 
This is located under `dependencies/scons/scons.py`, so simply run the following command:

.. code-block:: bash

  dependencies/scons/scons.py BUILD_TYPE=release

Because this is a long command, it is advisable to define a bash alias for this scons command. 
There are some predefined helper scripts that handle various frequently used compilation commands.
If you like, you can copy the following aliases to your `~/.bashrc` or `~/.bash_aliases` file, if you also have set the `OPENDIHU_HOME` environment variable as shown earlier.

.. code-block:: bash

  # define convenience commands for compilation
  alias scons='$OPENDIHU_HOME/dependencies/scons/scons.py'  
  alias s='scons'
  alias sd='$OPENDIHU_HOME/scripts/shortcuts/sd.sh'
  alias sdd='$OPENDIHU_HOME/scripts/shortcuts/sdd.sh'
  alias sddn='cd .. && scons BUILD_TYPE=d no_tests=yes no_examples=yes; cd -'
  alias sdn='scons BUILD_TYPE=d no_tests=yes no_examples=yes'
  alias srn='scons BUILD_TYPE=r no_tests=yes no_examples=yes'
  alias sr='$OPENDIHU_HOME/scripts/shortcuts/sr.sh'
  alias srd='$OPENDIHU_HOME/scripts/shortcuts/srd.sh'
  alias srr='$OPENDIHU_HOME/scripts/shortcuts/srr.sh'
  alias mkor='$OPENDIHU_HOME/scripts/shortcuts/mkor.sh'
  alias mkorn='$OPENDIHU_HOME/scripts/shortcuts/mkorn.sh'
  alias mkod='$OPENDIHU_HOME/scripts/shortcuts/mkod.sh'
  alias mkodn='$OPENDIHU_HOME/scripts/shortcuts/mkodn.sh'
  alias mkordn='$OPENDIHU_HOME/scripts/shortcuts/mkordn.sh'

Then, the following commands can be used for the build:

  * ``scons BUILD_TYPE=release`` or ``scons BUILD_TYPE=r`` or ``scons`` or ``s``:
    Build the file in the current directory in `release` mode, either to be used in the opendihu main directory to build the core library or in any example directory.
  
  * ``scons BUILD_TYPE=debug`` or ``scons BUILD_TYPE=d`` or ``sd``: Build `debug` target in current directory.
  * ``sdd``: To be used from within a `build_debug` directory. Go one directory up, build the example in `debug` target and go back to the original directory. This alias is equivalent to ``cd ..; scons BUILD_TYPE=debug; cd -``.
  * ``srr``: To be used from within a `build_release` directory. Go one directory up, build the example in `release` target and go back to the original directory. This alias is equivalent to ``cd ..; scons BUILD_TYPE=release; cd -``.
  * ``mkor``: "Make opendihu release". Use this command in any directory. It changes into the `opendihu` directory, executes `scons` there, to build the core library and changes back to the original directory.
  * ``mkorn``: "Make opendihu release, no tests". Same as `mkor`, except it does not build the unit tests. This is the most frequently used command to build the opendihu core.
  * ``mkod``: "Make opendihu debug". Use this command in any directory. It changes into the `opendihu` directory, executes `scons BUILD_TYPE=debug` there, to build the core library and changes back to the original directory.
  * ``mkodn``: "Make opendihu debug, no tests". Same as `mkor`, except it does not build the unit tests. This is the most frequently used command to build the opendihu core in debug target.
  * ``scons BUILD_TYPE=releasewithdebuginfo`` or ``scons BUILD_TYPE=rd`` or ``srd``: Build `releasewithdebuginfo` target in current directory.
  
As an example, if you work on a particular example and are in its `build_release` subdirectory, use ``mkorn && srr`` to build the core and the example and end up in the same directory afterwards.

If you have called `make` and everything has completed after some hours (green text), you were successful. Go on and build some examples (See next page, :doc:`getting_started`).
If not, read on, to find out what you need to configure in your case.

Configuring the build
^^^^^^^^^^^^^^^^^^^^^^^^

Configuration settings have to be provided in the python script `user-variables.scons.py`. These include settings for the dependency packages as well as further options concerning the build.

The option ``USE_VECTORIZED_FE_MATRIX_ASSEMBLY`` specifies if the Finite Element matrices should be assembled with SIMD instruction using the Vc library.
This leads to 4 elements always being a assembled at once using vector instructions (on systems with AVX-2).

If set to ``True``, this significantly speeds up the computation for problems that assemble a lot of matrices, e.g. solid mechanics problems.
However, it takes a long time to compile the code, up to 3x. If you intend to develop the core code, set it to ``False`` to have faster compilation. 
If you mainly want to run simulations including mechanics, set it to ``True``. 
(Also set it to ``False``, if compilation fails for ``True`` maybe because there is a bug somewhere that has not yet been found because the developers have this option always set to ``False``.)

For every dependency package there are variables like

.. code-block:: bash

  #PETSC_DOWNLOAD=True
  #PETSC_DIR="~/petsc/install"

(Note, `#` means commented out here, because you shouldn't specify both lines at once). 
The first line would instruct the build system to download and build the package, in this case PETSc. 
The second line would provide the path to an already existing installation on the system, which would then be used. Thus, specify either of those. 

There are similar options for all packages. You can read about more possibilities in the header of the `user-variables.scons.py` file. 

There are required dependencies, which need to be present in order for opendihu to work, and optional dependencies:

============================================================  ========  ===================================================================================
 Package                                                      Required    Description
============================================================  ========  ===================================================================================
`MPI`                                                             yes     | *Message Passing Interface*, used for data transfer between
                                                                          | processes. This should be your system MPI. If you let 
                                                                          | opendihu install it for you, `OpenMPI <https://www.open-mpi.org/>`_ 
                                                                          | will be chosen.
`PETSc <https://www.mcs.anl.gov/petsc/>`_                         yes     | Low-level data structures and solvers, see their `website <https://www.mcs.anl.gov/petsc/>`_
                                                                          | for more details.
`Python3`                                                         yes     | The `Python3 interpreter <https://www.python.org/>`_, 
                                                                          | version 3.9 or 3.6.5 for legacy. We need the development 
                                                                          | header and source files, therefore it is recommended to 
                                                                          | let opendihu build python for you, even if your system 
                                                                          | has python installed.
`pythonPackages`                                                  yes     | This is a custom collection of python packages for the
                                                                          | python 3 interpreter and are available in the
                                                                          | python configuration scripts. It consists of 
                                                                          | `numpy matplotlib scipy numpy-stl svg.path triangle geomdl vtk`.
`Base64 <https://github.com/tkislan/base64>`_                     yes     | An encoding standard and library that is used to create
                                                                          | binary VTK output files that can be viewed in Paraview.
                                                                          | Base64 encoded data is ASCII characters, the size is 4/3
                                                                          | of the raw binary data. The advantage is that despite 
                                                                          | being packed, it can be embedded in human-readable `XML`
                                                                          | files, which is the concept of VTK files.
`googletest <https://github.com/google/googletest>`_              no      | A testing framework, used for unit tests. Opendihu
                                                                          | compiles also without unit tests, but it is recommended 
                                                                          | to have them, especially for development of the core.
`SEMT <https://github.com/maierbn/semt>`_                         no      | This is a small C++ symbolic differentiation toolbox 
                                                                          | that will be used for nonlinear solid mechanics, to 
                                                                          | derive material laws.
`ADIOS2 <https://adios2.readthedocs.io/en/latest>`_               no      | Binary output file format and library, parallely 
                                                                          | efficient and self-descriptive. This is only installed, 
                                                                          | if you have a very recent version of `cmake`. If this
                                                                          | fails to install it is no problem as most users won't 
                                                                          | need it. It is needed for interfacing `MegaMol`.
`Vc <https://vcdevel.github.io/Vc-1.4.1/index.html>`_            yes      | A vectorization library that produces `simd` code 
                                                                          | depending on the hardware capabilities.
                                                                          |
`xbraid <https://github.com/XBraid/xbraid>`_                      no      | A framework for the parallel-in-time algorithm multigrid-
                                                                          | reduction-in-time (MGRIT)
`OpenCOR <https://opencor.ws/>`_                                  no      | `OpenCOR` is a modelling tool for CellML models and can 
                                                                          | convert `*.cellml` files to C code files, `*.c`. If
                                                                          | installed, the conversion of cellml input files is 
                                                                          | done automatically. If not, you can only input 
                                                                          | C files of the cellml models.
`libxml <http://xmlsoft.org/>`_                                    no     | A XML C parser, only needed for the installation of preCICE.
`preCICE <https://www.precice.org/>`_                              no     | Numerical coupling library, required, e.g., for the 
                                                                          | simulation of a muscle-tendon complex. This requires
                                                                          | a `boost <https://www.boost.org/>`_ installation as an additional prerequisite.
`Easylogging++ <https://github.com/zuhd-org/easyloggingpp>`_      yes     | This is the logging library. By default, logs are created 
                                                                          | in `/tmp/logs/` and output is written to the standard output.
============================================================  ========  =================================================================================== 

It is recommended to not let the build system download and build `MPI`, 
instead you should use your local MPI installation. 

On Ubuntu systems, the system MPI directory should already be set correctly by the default value in `user-variables.scons.py`. 
If you run `make`, you can check if MPI will be found.

If the MPI location is not detected automatically, you have to specify the path yourself. 
Find out in which path on your system MPI is installed. 
The required directory contains a `lib` and an `include` subdirectory. 
It may be located at `/usr/lib/openmpi`, `/usr/lib/mpich`, `/usr/lib/x86_64-linux-gnu/openmpi` or similar.
Set this path in `user-variables.scons.py` as the value of the variable `MPI_DIR`.

When running ``make``, ``make debug`` or ``make release``, the dependencies will be downloaded and installed, 
and consequently, debug or release target will be build. 
The installation of dependencies can take several hours. 
The compilation of the `core` afterwards completes in several minutes.

Troubleshooting
^^^^^^^^^^^^^^^^^^

If something fails during the installation, read the `config.log` file that will be created. 
It contains information about the commands used in the build process.

To restart the build process, it is sometimes required to clean the `scons` cache. This is done by deleting files ``.sconf_temp .sconsign.dblite`` which is executed by the command

.. code-block:: bash

  make clean

The dependencies that were already installed successfully will be detected the next time and not installed again. 
However, sometimes it is required to try to build a packages again.
You can force to rebuild selected packages by the `..._REBUILD` option, e.g.

.. code-block:: bash

  scons PETSC_REBUILD=True

to rebuild petsc, even if it was already detected. 

In general, the same options that can be specified in the `user-variables.scons.py` file 
can also be given like this on the command line as options to the `scons` command. (Also to the `sd` etc. shortcuts described earlier).

To restart with downloading the package and then installing it again, use the `..._REDOWNLOAD` option, like this:

.. code-block:: bash

  scons PETSC_REDOWNLOAD=True

Sometimes it also helps to delete the whole folder of a package in the `dependencies` subdirectory 
and retry the installation. 

If during execution of an example an error occurs that says numpy could not be imported, try to install the python packages of the python3 installation within opendihu yourself:

.. code-block:: bash

  opendihu/dependencies/python/install/bin/python3 -m pip install numpy matplotlib scipy svg.path geomdl

If a dependency fails to install, you can try to install it manually on your own. 
The commands that are used by the `scons` build system are printed to the console and additionally logged in the `config.log` file.

For advanced users, if you want to change the build system and update the commands that are executed
for installing a specific dependency, have a look at the directory `opendihu/dependencies/scons-config/sconsconfig/packages`.
It contains the source code for the build system. 
The main implementation is in `Package.py`, all other classes inherit from this class. 
Usually you find the file that is named like the dependency, e.g., `LAPACK.py` for Lapack or `PETSc.py` for PETSc.

If you change something here, you need to rebuild the python `egg` file of `scons-config`:

.. code-block:: bash

  cd <your-opendihu-path>
  cd dependencies/scons-config
  . install_manually.sh

Then, rerun the installation from the `opendihu` directory with `scons`.

If you don't succeed, ask for help and send us the `config.log` file.
