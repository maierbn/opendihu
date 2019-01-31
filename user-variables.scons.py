# Configuration for scons build system
#
# For each package the following variables are available:
# <PACKAGE>_DIR         Location of the package, must contain subfolders "include" and "lib" or "lib64" with header and library files.
# <PACKAGE>_INC_DIR     Location of (*.h) header files
# <PACKAGE>_LIB_DIR     Location of (*.a) libraries
# <PACKAGE>_LIBS        List of libraries, optional since the standard names are already hardcoded.
# <PACKAGE>_DOWNLOAD    Download, build and use a local copy of the package.
# <PACKAGE>_REDOWNLOAD  Force update of previously downloaded copy. For that <PACKAGE>_DOWNLOAD has to be also true.
# <PACKAGE>_REBUILD     Force a new build of the package without redownloading it if already has been downloaded earlier.
#
# You do one of the following:
# 1. Not specify any of the variables. Then standard locations in dependencies as well as /usr, /usr/local are searched.
# 2. Specify <PACKAGE>_DIR to directly give the base directory to the package's location. Do this to e.g. use system provided libraries.
# 3. Specify <PACKAGE>_INC_DIR and <PACKAGE>_LIB_DIR to point to the header and library directories. They are usually named "include" and "lib".
# 4. Set <PACKAGE>_DOWNLOAD=True or additionally <PACKAGE>_REDOWNLOAD=True to let the build system download and install everything on their own.

# set compiler to use
cc="gcc"   # c compiler
CC="g++"   # c++ compiler

# LAPACK, includes also BLAS, OpenBLAS is used
LAPACK_DOWNLOAD=True

# PETSc, this downloads and installs MUMPS (direct solver package) and its dependencies PT-Scotch, SCAlapack, ParMETIS, METIS
PETSC_DOWNLOAD=True

# Python 3.6
PYTHON_DOWNLOAD=True    # This downloads and uses Python, use it to be independent of an eventual system python

# Python packages - they are now all combined with the option PYTHONPACKAGES_DOWNLOAD
PYTHONPACKAGES_DOWNLOAD=True

# Base64, encoding library for binary vtk (paraview) output files
BASE64_DOWNLOAD=True

# Google Test, testing framework, not needed on Hazelhen
GOOGLETEST_DOWNLOAD=True

# SEMT, library for symbolic differentiation
SEMT_DOWNLOAD=True

# EasyLoggingPP, provides logging facilities
EASYLOGGINGPP_DOWNLOAD=True

# ADIOS2, adaptable I/O library, needed for interfacing MegaMol
ADIOS_DOWNLOAD=True

# MegaMol, visualization framework of VISUS, optional, needs ADIOS2
MEGAMOL_DOWNLOAD=False    # install MegaMol from official git repo, but needed is the private repo, ask for access to use MegaMol with opendihu

# MPI
# MPI is normally detected by runnig the mpicc command. If this is not available, you can provide the MPI_DIR as usual.
MPI_DIR="/usr/lib/openmpi"    # standard path for openmpi on ubuntu 16.04
#MPI_DIR="/usr/lib64/mpich/"

# automatically set MPI_DIR for other systems, like ubuntu 18.04 and Debian
try:
  import lsb_release
  lsb_info = lsb_release.get_lsb_information()   # get information about ubuntu version, if available
  if "RELEASE" in lsb_info:
    if lsb_info["RELEASE"] == "18.04":
      MPI_DIR="/usr/lib/x86_64-linux-gnu/openmpi"   # this is the standard path on ubuntu 18.04

  # use value of environment variable 'MPI_HOME' if it is set
  import os
  if os.environ.get("MPI_HOME") is not None:
    MPI_DIR = os.environ.get("MPI_HOME")
    
  # for Travis CI, build MPI ourselves
  if os.environ.get("TRAVIS") is not None:
    print "Travis CI detected, del MPI_DIR"
    del MPI_DIR
    MPI_DOWNLOAD=True
  
  # on neon use custom cmake
  if os.environ.get("HOSTNAME") == "neon":
    cmake="~/software/cmake/cmake-3.13.3-Linux-x86_64/bin/cmake"

except:
  pass

# download and build debugging MPI version
if False:
  del MPI_DIR
  MPI_DOWNLOAD=True
  MPI_IGNORE_MPICC=True    # this downloads and builds openmpi
  MPI_DEBUG=True            # this enables debugging flags such that valgrind memcheck can track MPI errors

# specialized settings for supercomputer (HazelHen)
import os
if os.environ.get("PE_ENV") is not None:  # if on hazelhen
  cc="cc"   # C compiler wrapper
  CC="CC"   # C++ compiler wrapper
  mpiCC="CC"  # mpi C++ compiler wrapper

  # use cray-pat for profiling
  USE_CRAY_PAT=False

  # use -hpl option with cray compiler to create an optimization program library
  USE_HPL=False

  # do not use googletest
  GOOGLETEST_DOWNLOAD=False  

  # do not use buggy python packages
  PYTHONPACKAGES_DOWNLOAD=False

  #MPI_DIR = os.environ.get("CRAY_MPICH_DIR")
  #LAPACK_DOWNLOAD = False
  #LAPACK_DIR = os.environ.get("CRAY_LIBSCI_PREFIX_DIR")
  #PETSC_DOWNLOAD = False
  #PETSC_DIR = os.environ.get("PETSC_DIR")

# Steps for getting started on HazelHen
#   module swap PrgEnv-cray/6.0.4 PrgEnv-gnu  # to switch to GNU programming environment, however also Intel and Cray environments work
#   module load cray-libsci
#   module load cray-petsc  (or cray-petsc-64 for big data)




