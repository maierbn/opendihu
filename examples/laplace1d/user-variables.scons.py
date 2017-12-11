# Configuration for scons build system
#
# For each package the following variables are available:
# <PACKAGE>_DIR         Location of the package, must contain subfolders "include" and "lib" or "lib64" with header and library files.
# <PACKAGE>_INC_DIR     Location of (*.h) header files
# <PACKAGE>_LIB_DIR     Location of (*.a) libraries
# <PACKAGE>_LIBS        List of libraries, optional since the standard names are already hardcoded.
# <PACKAGE>_DOWNLOAD    Download, build and use a local copy of the package.
# <PACKAGE>_REDOWNLOAD  Force update of previously downloaded copy. For that <PACKAGE>_DOWNLOAD has to be also true.
#
# You do one of the following:
# 1. Not specify any of the variables. Then standard locations in dependencies as well as /usr, /usr/local are searched
# 2. Specify <PACKAGE>_DIR to directly give the base directory to the package's location
# 3. Specify <PACKAGE>_INC_DIR and <PACKAGE>_LIB_DIR to point to the header and library files
# 4. Set <PACKAGE>_DOWNLOAD=True or additionally <PACKAGE>_REDOWNLOAD=True to let the build system download and install everything on their own.


# MPI
#MPI_DIR="/usr/lib/openmpi"

# LAPACK, includes also BLAS
#LAPACK_DOWNLOAD=False
#LAPACK_REDOWNLOAD=False      # force download, even if the files are already present

# PETSc
#PETSC_DOWNLOAD=False
#PETSC_REDOWNLOAD=False
#PETSC_DIR="../dependencies/petsc"
