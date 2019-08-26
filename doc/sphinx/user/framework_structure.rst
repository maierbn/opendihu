Framework structure
=======================

In the following a first overview of the structure of opendihu is given.

Usage
------

To use the framework, one needs to use two programming languages: C++ and Python.
For a simulation, one has to write a small C++ `main` file that includes the main framework header (`opendihu.h`) and contains a definition of the desired nested solvers as class templates.

This source file will be compiled either in debug or release mode. It will be linked against a single opendihu library under `core/build_debug/libopendihud.a` or `core/build_debug/libopendihu.a`.
For example, you can build any example in its top-level directory by calling 

.. code-block:: bash

  scons BUILD_TYPE=d

The resulting executable in this case will be created under the `build_debug` subdirectory. It can be executed with the file name of a python script as its argument. This python script contains all the settings.
For example, if the created executable is `build_debug/simulation` and the python settings file is `settings.py`, one would typically execute the simulation in the following way:

.. code-block:: bash

  cd build_debug
  ./simulation ../settings.py

Depending on the settings, e.g., a new folder `out` will be created with lots of output files.

To change the settings, you have to edit the python settings file and rerun the simulation. There is no need to recompile everything.

Usually, large simulations are run in parallel. Most simulations can directly run with multiple cores, without any changes in the settings. To execute the example with 4 processes, simply run

.. code-block:: bash

  cd build_debug
  mpirun -n 4 ./simulation ../settings.py

To get more information, read :doc:`/user/getting_started`.


Directory structure
---------------------

The `opendihu` contains the following subdirectories:

=================   =================================
core                 | This contains the C++ code of the framework, implementing all functionality. 
                     | Users usually don't have to touch this directory, developers will.
                     | The `core/src` contains the actual code, the binary library files will be build
                     | under `core/build_debug` and `core/build_release`.
dependencies         | Here all external dependencies are stored and installed in an own folder for each.
                     | This directory will be populated by the build system `scons` that can download
                     | and build all dependencies.
doc                  | This directory is for documentation and contains useful pdfs, some own theory 
                     | documentation under `doc/derivations` and some symbolic computations using the
                     | python symbolic toolbox `sympy`.
examples              This contains all examples that use the framework. Users should work in this directory.
scripts              | This is a collection of useful Python scripts, such as the plot script or a
                     | reader utility for the python output format of opendihu.
testing              This contains unit tests and system tests.
=================   =================================

