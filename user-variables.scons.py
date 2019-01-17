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

# set compiler
cc="gcc"   # c compiler
CC="g++"   # c++ compiler

# LAPACK, includes also BLAS, current OpenBLAS is used
LAPACK_DOWNLOAD=True

# PETSc
PETSC_DOWNLOAD=True

# Python
PYTHON_DOWNLOAD=True    # This downloads and uses Python, use it to be independent of an eventual system python

# Python packages - they are now all combined with the option PYTHONPACKAGES_DOWNLOAD
# Numpy
#CYTHON_DOWNLOAD=True
#NUMPYC_DOWNLOAD=True

# SciPy
#SCIPY_DOWNLOAD=True

# Matplotlib and other python dependencies
#BZIP2_DOWNLOAD=True
#MATPLOTLIB_DOWNLOAD=False
#NUMPYSTL_DOWNLOAD=False
#SVGPATH_DOWNLOAD=False

PYTHONPACKAGES_DOWNLOAD=True

# Base64
BASE64_DOWNLOAD=True

# Google Test
GOOGLETEST_DOWNLOAD=True

# SEMT
SEMT_DOWNLOAD=True

# EasyLoggingPP
EASYLOGGINGPP_DOWNLOAD=True

# MegaMol
#MEGAMOL_DOWNLOAD=True     # install MegaMol from official git repo
#USE_MEGAMOL=True          # link to mmconsole main function and start megamol, if config["MegaMolArguments"] is set

# MPI
# MPI is normally detected using mpicc. If this is not available, you can provide the MPI_DIR as usual.
MPI_DIR="/usr/lib/openmpi"    # standard path for ubuntu 16.04
#MPI_DIR="/usr/lib64/mpich/"

# automatically set MPI_DIR for ubuntu 18.04
try:
  import lsb_release
  lsb_info = lsb_release.get_lsb_information()   # get information about ubuntu version, if available
  if "RELEASE" in lsb_info:
    if lsb_info["RELEASE"] == "18.04":
      MPI_DIR="/usr/lib/x86_64-linux-gnu/openmpi"   # this is the path for ubuntu 18.04

  # use value of environment variable 'MPI_HOME' if it is set
  import os
  if os.environ.get("MPI_HOME") is not None:
    MPI_DIR = os.environ.get("MPI_HOME")
    
  # for Travis CI build MPI ourselves
  if os.environ.get("TRAVIS") is not None:
    print "Travis CI detected, del MPI_DIR"
    del MPI_DIR
    MPI_DOWNLOAD=True

except:
  pass

# download and build debugging MPI version
if False:
  del MPI_DIR
  MPI_DOWNLOAD=True
  MPI_IGNORE_MPICC=True    # this downloads and builds openmpi
  MPI_DEBUG=True            # this enables debugging flags such that valgrind memcheck can track MPI errors

# other variables for hazelhen
import os
if os.environ.get("PE_ENV") is not None:  # if on hazelhen
  cc="cc"   # C compiler
  CC="CC"   # C++ compiler
  mpiCC="CC"  # mpi C++ compiler

  # use cray-pat for profiling
  USE_CRAY_PAT=False

  # use -hpl option with cray compiler to create an optimization program library
  USE_HPL=False

  #MPI_DIR = os.environ.get("CRAY_MPICH_DIR")
  #LAPACK_DOWNLOAD = False
  #LAPACK_DIR = os.environ.get("CRAY_LIBSCI_PREFIX_DIR")
  #PETSC_DOWNLOAD = False
  #PETSC_DIR = os.environ.get("PETSC_DIR")

# module restore opendihu
# or 
#   module swap PrgEnv-cray/6.0.4 PrgEnv-gnu
#   module load cray-libsci
#   module load cray-petsc  (or cray-petsc-64 for big data)




