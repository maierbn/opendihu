import sys, os
import socket
from distutils import sysconfig
from Package import Package

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

def parse_conf(ctx, conf_path, lib_dirs, libs):
  vars = {}
  sysconfig.parse_makefile(conf_path, vars)
  flag_dict = ctx.env.ParseFlags(vars['PACKAGES_LIBS'])
  lib_dirs.extend(flag_dict['LIBPATH'])
  for ii in range(len(libs)):
    libs[ii].extend(flag_dict['LIBS'])

def find_conf(ctx, base, inc_dirs, lib_dirs, libs, extra_libs):
  # PETSc 3.1
  conf_path = os.path.join(base, 'conf', 'petscvariables')
  if os.path.exists(conf_path):
      parse_conf(ctx, conf_path, lib_dirs, libs)

  # PETSC 2.3.3
  conf_path = os.path.join(base, 'bmake', 'petscconf')
  if os.path.exists(conf_path):
    vars = {}
    sysconfig.parse_makefile(conf_path, vars)
    if 'PETSC_ARCH' in vars:
        arch = vars['PETSC_ARCH']
        inc_dirs.extend([os.path.join(base, 'bmake', arch)])
        lib_dirs.extend([os.path.join(base, 'lib', arch)])
        conf_path = os.path.join(base, 'bmake', arch, 'petscconf')
        parse_conf(ctx, conf_path, lib_dirs, libs)

class PETSc(Package):

  def __init__(self, **kwargs):
    defaults = {
      #'download_url': 'http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.6.tar.gz',
      'download_url': 'http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.12.3.tar.gz',
    }
    defaults.update(kwargs)
    super(PETSc, self).__init__(**defaults)
    self.sub_dirs = [('include','lib')]
    self.libs = [['cmumps', 'HYPRE', 'sundials_cvode', 'petsc', 'scalapack', 'parmetis', 'dmumps', 'smumps', 'zmumps', 'mumps_common', 'cmumps', 'scalapack', 'petsc', 'pord', 'parmetis', 'petsc', 'dmumps', 'sundials_nvecparallel', 'sundials_nvecserial', 'petsc', 'sundials_cvode', 'smumps']]
      # ['petsc', 'cmumps', 'dmumps', 'HYPRE', 'mumps_common', 'pord', 'scalapack', 'smumps', 'sundials_cvode', 'sundials_nvecparallel', 'sundials_nvecserial', 'zmumps', 'parmetis']]
    self.headers = ['petsc.h']

    self.check_text = petsc_text
    self.static = False
    
    if os.environ.get("PE_ENV") is not None:  # if on hazelhen
      print("Same for Petsc.")
    
      # on hazel hen login node do not run MPI test program because this is not possible (only compile)
      self.run = False
      
    self.number_output_lines = 4121
      
  def check(self, ctx):
    if os.environ.get("PE_ENV") is not None:  # if on hazelhen
      ctx.Message('Not checking for PETSc ... ')
      ctx.Result(True)
      return True
  
    env = ctx.env
    
    # --with-cc='+env["CC"]+'\
    
    # debugging build handler 
    if self.have_option(env, "PETSC_DEBUG"):
      # debug build with MUMPS
      print("PETSc debugging build is on!")
      self.set_build_handler([
        'mkdir -p ${PREFIX}',
        './configure --prefix=${PREFIX} --with-debugging=yes --with-shared-libraries=1 \
        --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
          ---with-cc='+env["mpicc"]+'\
        --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch --download-sundials --download-hypre \
         | tee out.txt',
        '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
        '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        'ln -fs ${PREFIX}/lib/libparmetis.so ${PREFIX}/lib/parmetis.so'    # create parmetis.so link for chaste
      ])
    else:
      # standard release build with MUMPS
      # This needs bison installed
      
      # for metis to work, we need --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch        
      self.set_build_handler([
          'mkdir -p ${PREFIX}',
          #'PATH=${PATH}:${DEPENDENCIES_DIR}/bison/install/bin \
          './configure --prefix=${PREFIX} --with-debugging=no --with-shared-libraries=1 \
          --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
          ---with-cc='+env["mpicc"]+'\
          --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch --download-sundials --download-hypre \
          COPTFLAGS=-O3\
          CXXOPTFLAGS=-O3\
          FOPTFLAGS=-O3 | tee out.txt',
         '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',     # do it twice, the first time fails with PGI
         '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
         '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
         '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
         'ln -fs ${PREFIX}/lib/libparmetis.so ${PREFIX}/lib/parmetis.so'    # create parmetis.so link for chaste
      ])
    
    ctx.Message('----------------------------------------------------\nNote that PETSc has been updated to version 3.12.3. \nTo update, run \'scons PETSC_REDOWNLOAD=True\'.\n(This message is independent of the currently installed version.)\n----------------------------------------------------\n')
    ctx.Message('Checking for PETSc ...         ')
    self.check_options(env)

    res = super(PETSc, self).check(ctx, loc_callback=find_conf)
    #self.check_required(res[0], ctx)
  
    # if installation of petsc fails, retry without mumps and extra packages like parmetis, hdf5 or hypre
    if not res[0] and socket.gethostname()!= 'cmcs09':
      ctx.Log('Retry without MUMPS\n')
      ctx.Message('Retry to install a fail-back PETSc without MUMPS, Hypre, SUNDIALS and ParMETIS ...')
      if "PETSC_REDOWNLOAD" in Package.one_shot_options:
        Package.one_shot_options.remove('PETSC_REDOWNLOAD')
      if "PETSC_REBUILD" in Package.one_shot_options:
        Package.one_shot_options.remove('PETSC_REBUILD')
      
      if self.have_option(env, "PETSC_DEBUG"):
        # debug build, without MUMPS
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=yes \
            --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
            --with-cc='+env["mpicc"]+' | tee out.txt',
          '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
          '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
        ])
      else:
        # release build without MUMPS
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
          --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
          --with-cc='+env["mpicc"]+'\
          COPTFLAGS=-O3\
          CXXOPTFLAGS=-O3\
          FOPTFLAGS=-O3 | tee out.txt',
        '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt || make',
        '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt) || make install',
        ])
      self.libs = ['petsc']
      
      self.number_output_lines = 3990
      
      res = super(PETSc, self).check(ctx, loc_callback=find_conf)
      self.check_required(res[0], ctx)
    
    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
