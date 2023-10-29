import sys, os
import socket
from .Package import Package

petsc_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <petsc.h>
int main(int argc, char* argv[]) {
   PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
   printf("MPI version %d.%d\n", MPI_VERSION, MPI_SUBVERSION);
   printf("Petsc version %d.%d.%d\n", PETSC_VERSION_MAJOR,PETSC_VERSION_MINOR,PETSC_VERSION_SUBMINOR);
   PetscFinalize();
   return EXIT_SUCCESS;
}
'''

class PETSc(Package):

  def __init__(self, **kwargs):
    defaults = {
      'download_url': 'http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.15.5.tar.gz'
    }
    defaults.update(kwargs)
    super(PETSc, self).__init__(**defaults)
    self.sub_dirs = [('include','lib')]
    self.libs = [['cmumps', 'HYPRE', 'petsc', 'scalapack', 'parmetis', 'dmumps', 'smumps', 'zmumps', 'mumps_common', 'cmumps', 'scalapack', 'petsc', 'pord', 'parmetis', 'petsc', 'dmumps', 'petsc', 'smumps']]
      # ['petsc', 'cmumps', 'dmumps', 'HYPRE', 'mumps_common', 'pord', 'scalapack', 'smumps', 'sundials_cvode', 'sundials_nvecparallel', 'sundials_nvecserial', 'zmumps', 'parmetis']]
    self.headers = ['petsc.h']

    self.check_text = petsc_text
    self.static = False
    
    if os.environ.get("PE_ENV") is not None:  # if on hazelhen
      print("Same for Petsc.")
    
      # on hazel hen login node do not run MPI test program because this is not possible (only compile)
      self.run = False
      
    self.number_output_lines = 4360
      
  def check(self, ctx):
    if os.environ.get("PE_ENV") is not None:  # if on hazelhen
      ctx.Message('Not checking for PETSc ... ')
      ctx.Result(True)
      return True
  
    env = ctx.env
    ctx.Message('Checking for PETSc ...         ')
    
    # --with-cc='+env["CC"]+'\
    # --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
    # --with-cc='+env["mpicc"]+'\


    # special case for host cmcs05
    # nur für development der gpu-isierung benötigt (üblicherweise auf Rechner cmcs05): 
    if socket.gethostname() == "cmcs05":
      self.set_build_handler([
        'mkdir -p ${PREFIX}',
        '$rm -rf $(ls | grep "linux")',
        #'PATH=${PATH}:${DEPENDENCIES_DIR}/bison/install/bin \
        './configure --prefix=${PREFIX} --with-debugging=no --with-shared-libraries=1 \
          --download-fblaslapack=1 \
          --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch --download-sundials --download-hypre \
          COPTFLAGS=-O3\
          CXXOPTFLAGS=-O3\
          --with-mpi-dir=${MPI_DIR} --with-batch\
          FOPTFLAGS=-O3 | tee out.txt',
        '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',     # do it twice, the first time fails with PGI
        '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
        '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
      ])
      res = super(PETSc, self).check(ctx)
      self.check_required(res[0], ctx)
      ctx.Result(res[0])
      return res[0]


    if self.have_option(env, "PETSC_DEBUG"):
      # standard debug build, without any mpi related option
      print("PETSc debugging build is on!")
      self.set_build_handler([
        'mkdir -p ${PREFIX}',
        '$rm -rf $(ls | grep "linux")',
        './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=yes \
          --download-fblaslapack=1 \
          --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch --download-sundials --download-hypre \
         | tee out.txt',
        '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
        '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        'ln -fs ${PREFIX}/lib/libparmetis.so ${PREFIX}/lib/parmetis.so'    # create parmetis.so link for chaste
      ]) #  --with-batch benötig
    else:
      # standard release build with MUMPS, without any mpi related option
      # This needs bison installed
      
      # for metis to work, we need --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch
      self.set_build_handler([
        'mkdir -p ${PREFIX}',
        '$rm -rf $(ls | grep "linux")',
        './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no  \
          --download-fblaslapack=1 \
          --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch --download-sundials --download-hypre \
          COPTFLAGS=-O3\
          CXXOPTFLAGS=-O3\
          FOPTFLAGS="-O3 -fallow-argument-mismatch" | tee out.txt',
        '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',     # do it twice, the first time fails with PGI
        '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
        '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        'ln -fs ${PREFIX}/lib/libparmetis.so ${PREFIX}/lib/parmetis.so'    # create parmetis.so link for chaste
      ])
    
    self.check_options(env)
    res = super(PETSc, self).check(ctx)
  

    # if installation of petsc failed with the current command, retry with different options 
    if not res[0]:
      ctx.Log('Retry (1) with option --with-cc=\n')
      ctx.Message('Retry (1) using --with-cc ...')
      if "PETSC_REDOWNLOAD" in Package.one_shot_options:
        Package.one_shot_options.remove('PETSC_REDOWNLOAD')
      if "PETSC_REBUILD" in Package.one_shot_options:
        Package.one_shot_options.remove('PETSC_REBUILD')
      
      if self.have_option(env, "PETSC_DEBUG"):
        # debug build, without --with-mpi-dir, but with --with-cc
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          '$rm -rf $(ls | grep "linux")',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=yes \
            --with-cc='+env["mpicc"]+' --with-cxx='+env["mpiCC"]+' --with-batch \
            --download-fblaslapack=1 \
            --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch --download-sundials --download-hypre \
           | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        ])
      else:
        # release build, without --with-mpi-dir, but with --with-cc
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          '$rm -rf $(ls | grep "linux")',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
            --with-cc='+env["mpicc"]+' --with-cxx='+env["mpiCC"]+' \
            --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch --download-sundials --download-hypre \
            COPTFLAGS=-O3 \
            CXXOPTFLAGS=-O3 \
            FOPTFLAGS="-O3 -fallow-argument-mismatch" | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt || make',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt) || make install',
        ])  

      # check again, if this works
      self.check_options(env)
      res = super(PETSc, self).check(ctx)


    # if it also did not work
    if not res[0]:
      ctx.Log('Retry (2) with option --with-mpi-dir\n')
      ctx.Message('Retry (2) with option --with-mpi-dir ...')
      if self.have_option(env, "PETSC_DEBUG"):
        # debug build, with --with-mpi-dir, without --with-cc
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          '$rm -rf $(ls | grep "linux")',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=yes \
            --with-mpi-dir=${MPI_DIR} --with-batch \
            --download-fblaslapack=1 \
            --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch --download-sundials --download-hypre \
           | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        ])
      else:
        # release build, with --with-mpi-dir, without --with-cc
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          '$rm -rf $(ls | grep "linux")',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
            --with-mpi-dir=${MPI_DIR} \
            --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch --download-sundials --download-hypre \
            COPTFLAGS=-O3 \
            CXXOPTFLAGS=-O3 \
            FOPTFLAGS="-O3 -fallow-argument-mismatch" | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt || make',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt) || make install',
        ])  

      # try again
      self.check_options(env)
      res = super(PETSc, self).check(ctx)


    # if it also did not work
    if not res[0]: 
      # if installation of petsc failed again, retry without mumps and extra packages like parmetis, hdf5 or hypre
      ctx.Log('Retry (3) without MUMPS, Hypre, SUNDIALS and ParMETIS\n')
      ctx.Message('Retry (3) without MUMPS, Hypre, SUNDIALS and ParMETIS ...')
      if "PETSC_REDOWNLOAD" in Package.one_shot_options:
        Package.one_shot_options.remove('PETSC_REDOWNLOAD')
      if "PETSC_REBUILD" in Package.one_shot_options:
        Package.one_shot_options.remove('PETSC_REBUILD')
        
      if self.have_option(env, "PETSC_DEBUG"):
        # debug build, without MUMPS
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          '$rm -rf $(ls | grep "linux")',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=yes \
            --with-batch \
            --download-fblaslapack=1 | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        ])
      else:
        # release build without MUMPS
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          '$rm -rf $(ls | grep "linux")',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
            COPTFLAGS=-O3 \
            CXXOPTFLAGS=-O3 \
            FOPTFLAGS=-O3 | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt || make',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt) || make install',
        ])
      self.libs = ['petsc']
      self.number_output_lines = 3990
    
      # try again
      self.check_options(env)
      res = super(PETSc, self).check(ctx)


    # if it also did not work
    if not res[0]:
      ctx.Log('Retry (4) without MUMPS, Hypre, SUNDIALS and ParMETIS, with --with_mpi-dir\n')
      ctx.Message('Retry (4) without MUMPS, with --with_mpi-dir ...')
      if "PETSC_REDOWNLOAD" in Package.one_shot_options:
        Package.one_shot_options.remove('PETSC_REDOWNLOAD')
      if "PETSC_REBUILD" in Package.one_shot_options:
        Package.one_shot_options.remove('PETSC_REBUILD')

      if self.have_option(env, "PETSC_DEBUG"):
        # debug build, without MUMPS
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          '$rm -rf $(ls | grep "linux")',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=yes \
            --with-mpi-dir=${MPI_DIR} --with-batch \
            --download-fblaslapack=1 | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        ])
      else:
        # release build without MUMPS
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          '$rm -rf $(ls | grep "linux")',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
            --with-mpi-dir=${MPI_DIR} \
            COPTFLAGS=-O3 \
            CXXOPTFLAGS=-O3 \
            FOPTFLAGS=-O3 | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt || make',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt) || make install',
        ])

      self.check_options(env)
      res = super(PETSc, self).check(ctx)
  

    # if it also did not work
    if not res[0]:
      ctx.Log('Retry (5) without MUMPS, Hypre, SUNDIALS and ParMETIS, Fortran bindings, with --with_mpi-dir, f2cblaslapack (also in release)\n')
      ctx.Message('Retry (5) without MUMPS, Hypre, SUNDIALS and ParMETIS, Fortran bindings, with --with_mpi-dir, f2cblaslapack (also in release) ...')
      # --download-fblaslapack needs Fortran, use --download-f2cblaslapack instead
      if "PETSC_REDOWNLOAD" in Package.one_shot_options:
        Package.one_shot_options.remove('PETSC_REDOWNLOAD')
      if "PETSC_REBUILD" in Package.one_shot_options:
        Package.one_shot_options.remove('PETSC_REBUILD')

      if self.have_option(env, "PETSC_DEBUG"):
        # debug build, without MUMPS
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          '$rm -rf $(ls | grep "linux")',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=yes \
            --with-mpi-dir=${MPI_DIR} --with-batch \
            --with-fortran-bindings=0 --with-fc=0 \
            --download-f2cblaslapack=1 | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        ])
      else:
        # release build without MUMPS
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          '$rm -rf $(ls | grep "linux")',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
            --with-mpi-dir=${MPI_DIR} \
            --with-fortran-bindings=0 --with-fc=0 \
            --download-f2cblaslapack=1 \
            COPTFLAGS=-O3 \
            CXXOPTFLAGS=-O3 \
            FOPTFLAGS=-O3 | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt || make',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt) || make install',
        ])

      self.check_options(env)
      res = super(PETSc, self).check(ctx)
    

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
