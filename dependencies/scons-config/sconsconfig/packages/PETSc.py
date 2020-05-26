import sys, os
from distutils import sysconfig
from Package import Package

petsc_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <petsc.h>
int main(int argc, char* argv[]) {
   PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
   printf("%d\n", MPI_VERSION);
   printf("%d\n", MPI_SUBVERSION);
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
            'download_url': 'http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.13.1.tar.gz',
        } #http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.6.tar.gz
        defaults.update(kwargs)
        super(PETSc, self).__init__(**defaults)
        self.sub_dirs = [('include','lib')]
        self.libs = [['petsc'], ['petscksp', 'petscvec', 'petsc']]

        if os.environ.get("SITE_PLATFORM_NAME") == "hazelhen":
          if os.environ.get("PE_ENV") == "GNU":
            self.libs = ["craypetsc_gnu_real"]
            self.extra_libs = ["sci_gnu_71_mpi_mp"]
            print("{} environment detected, using \"{}\" for Petsc".format(os.environ.get("PE_ENV"), self.libs[0]))
          else:
            print("WARNING: The PE environment seems to be {}, not GNU, this is not supported".format(os.environ.get("PE_ENV")))
        
          # on hazel hen login node do not run MPI test program because this is not possible (only compile)
          self.run = False
          
        self.check_text = petsc_text
        self.static = False
        
        # Setup the build handler. This needs bison installed.

# # # # # # comment, muss eingearbeitet werden, so dass automatisch benutzt wird, wenn pgi benutzt wird:
# //cmcs05 'ist': --with-mpi-dir...    ...=/usr/lib/x86_64-linux-gnu/openmpi
# //      'soll': ...=/afs/mathematik.uni-stuttgart.de/home/kraemer/share/environment-modules/Packages/openmpi/3.0.0/gcc4.9.x
# ./configure --prefix=/usr/local/home/kraemer/software/opendihu2/dependencies/petsc/install --with-shared-libraries=1 --with-debugging=no --with-blas-lapack-lib=/usr/local/home/kraemer/software/opendihu2/dependencies/lapack/install/lib/libopenblas.so --with-mpi-dir=/afs/mathematik.uni-stuttgart.de/home/kraemer/share/environment-modules/Packages/openmpi/3.0.0/gcc4.9.x --download-mumps --download-scalapack --download-parmetis --download-metis COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3 --with-batch CC=$CC CXX=$CXX FC=$FC

        self.set_build_handler([
            #'PATH=${PATH}:${DEPENDENCIES_DIR}/bison/install/bin \
            './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
            --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
            --with-mpi-dir=${MPI_DIR}\
            --download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch \
            COPTFLAGS=-O3\
            CXXOPTFLAGS=-O3\
            FOPTFLAGS=-O3\
            --with-batch',  # FC, CC und CXX nicht angeben, er nimmt sie sowieso, sagt nur er w+rde ignorieren, hat das aber schon über nen anderen Kanal....
            'make all',     # do not add -j option, because it is not supported by Makefile of PETSc
            'echo "sleep 3 s" && sleep 3',
            'make install',
            'make test',
        ])
        
        self.number_output_lines = 4121
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for PETSc ... ')
        self.check_options(env)

        res = super(PETSc, self).check(ctx, loc_callback=find_conf)
        
        # if installation of petsc fails, retry without mumps
        if not res[0]:
          ctx.Log('Retry without MUMPS')
          ctx.Message('Retry to install PETSc without MUMPS ...')
          
          # Setup the build handler.
          self.set_build_handler([
              './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
              --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
              --with-mpi-dir=${MPI_DIR}\
              COPTFLAGS=-O3\
              CXXOPTFLAGS=-O3\
              FOPTFLAGS=-O3',
              'make all',     # do not add -j option, because it is not supported by Makefile of PETSc
              'echo "sleep 3 s" && sleep 3',
              'make install',
              'make test',
          ])
          
          self.number_output_lines = 3990
          
          self.check_options(env)

          res = super(PETSc, self).check(ctx, loc_callback=find_conf)

          self.check_required(res[0], ctx)
          ctx.Result(res[0])
          
        return res[0]
