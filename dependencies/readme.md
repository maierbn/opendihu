This directory contains all external dependency packages. The subdirectory `scons` contains the scons build system which is written in python. The subdirectory `scons-config` contains further python classes used to download and install dependencies. All other subdirectories correspond to an external dependency each. They can be deleted because the build system is able to download and build all dependencies.

The directory structure for a dependency is as follows:
dependencyname - lower case name of the package
  src          - contains all extracted files of the package
  install      - contains the installed interface that can be used by the core 
    lib        - contains libraries, this is added to the library path of the linker
    include    - contains header files to be included, this is added to the include path of the compiler
