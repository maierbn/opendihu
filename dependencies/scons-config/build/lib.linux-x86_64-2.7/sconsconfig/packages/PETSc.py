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
            'download_url': 'http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.6.tar.gz',
        }
        defaults.update(kwargs)
        super(PETSc, self).__init__(**defaults)
        self.sub_dirs = [('include','lib')]
        self.libs = [['petsc', 'cmumps', 'dmumps', 'HYPRE', 'mumps_common', 'pord', 'ptesmumps', 'scalapack', 'scotch', 'scotcherr',
          'scotcherrexit', 'smumps', 'sundials_cvode', 'sundials_nvecparallel', 'sundials_nvecserial', 'zmumps']]
          

        self.check_text = petsc_text
        self.static = False
        
        if os.environ.get("PE_ENV") is not None:  # if on hazelhen
          
          #if os.environ.get("PE_ENV") == "GNU":
          #  self.libs = ["craypetsc_gnu_real"]
          #  self.extra_libs = ["sci_gnu_71_mpi_mp"]
          #  print("{} environment detected, using \"{}\" for Petsc".format(os.environ.get("PE_ENV"), self.libs[0]))
          #else:
          #  print("WARNING: The PE environment seems to be {}, not GNU, this is not supported".format(os.environ.get("PE_ENV")))
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
        
        
        # debugging build handler 
        if self.have_option(env, "PETSC_DEBUG"):
          # debug build with MUMPS
          print("PETSc debugging build is on!")
          self.set_build_handler([
            'mkdir -p ${PREFIX}',
            #'PATH=${PATH}:${DEPENDENCIES_DIR}/bison/install/bin \
            './configure --prefix=${PREFIX} --with-debugging=yes \
            --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
            --with-mpi-dir=${MPI_DIR}\
            --download-mumps --download-scalapack --download-parmetis --download-metis | tee out.txt',
            '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
            '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
            'ln -fs ${PREFIX}/lib/libparmetis.so ${PREFIX}/lib/parmetis.so'    # create parmetis.so link for chaste
          ])
        else:
          # standard release build with MUMPS
          # This needs bison installed
          if socket.gethostname() != 'cmcs09':
            # on normal host
            
            self.set_build_handler([
                'mkdir -p ${PREFIX}',
                #'PATH=${PATH}:${DEPENDENCIES_DIR}/bison/install/bin \
                './configure --prefix=${PREFIX} --with-debugging=no \
                --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
                --with-mpi-dir=${MPI_DIR}\
                --download-mumps --download-scalapack --download-parmetis --download-metis --download-sundials --download-hypre \
                COPTFLAGS=-O3\
                CXXOPTFLAGS=-O3\
                FOPTFLAGS=-O3 | tee out.txt',
               '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
               '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
               'ln -fs ${PREFIX}/lib/libparmetis.so ${PREFIX}/lib/parmetis.so'    # create parmetis.so link for chaste
            ])
          else:                              
            # on cmcs09 using PGI
                                                             # # # # # P G I # # # # #
            #print("WARNING: MPI_DIR is set manually in scons-config/sconsconfig/packages/PETSc.Py." ) # because --with-mpi-dir=${MPI_DIR} does not work
            print("INFO: setting FLAG '--with-mpiexec' manually in PETSc.Py. ")
            self.set_build_handler([ 
                'mkdir -p ${PREFIX}',
 # don't use CC=$CC nor CXX=$CXX such that compiler can choose mpicc and mpicxx instead
 # --with-mpi=0 -I/usr/local/home/kraemer/opendihu/dependencies/petsc/install/include/petsc/mpiuni\
 # --with-mpi-include=/usr/local/home/kraemer/offloading/pgi_gcc7.2.0/linux86-64/2018/mpi/openmpi-2.1.2/include \ can't use both include and dir.
 # ---with-mpiexec=/usr/local/home/kraemer/offloading/pgi_gcc7.2.0/linux86-64/2018/mpi/openmpi/bin/mpirun\               
 #'PATH=${PATH}:${DEPENDENCIES_DIR}/bison/install/bin \
 #--CCFLAGS="-I/usr/local/home/kraemer/offloading/pgi_gcc7.2.0/linux86-64/2018/mpi/openmpi-2.1.2/include" \
 #--CFLAGS="-I/usr/local/home/kraemer/offloading/pgi_gcc7.2.0/linux86-64/2018/mpi/openmpi-2.1.2/include" \
 #--CFLAGS="-L/afs/.mathe/home/cmcs/share/environment-modules/Packages/gcc/7.2.0/lib/gcc/x86_64-pc-linux-gnu/7.2.0"\ #might be needed otherwise gcc4.9 libs might end up in config
                './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
                --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so \
                --with-mpi-dir=${MPI_DIR} \
                --with-fc=0 \
                --with-mpiexec=/usr/local/home/kraemer/offloading/pgi_gcc7.2.0/linux86-64/2018/mpi/openmpi/bin/mpirun \
                COPTFLAGS=-fast \
                CXXOPTFLAGS=-fast | tee out.txt',
               '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
               '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
               #'cp /usr/local/home/kraemer/offloading/pgi_gcc7.2.0/linux86-64/2018/mpi/openmpi-2.1.2/include/mpi.h /usr/local/home/kraemer/opendihu/dependencies/petsc/install/include/',
               #'ln -sfn /usr/local/home/kraemer/offloading/pgi_gcc7.2.0/linux86-64/2018/mpi/openmpi-2.1.2/include /usr/local/home/kraemer/opendihu/dependencies/petsc/install/include/mpiinclude',
            ])
        
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
          
          # Setup the build handler.
          
          
          if self.have_option(env, "PETSC_DEBUG"):
            # debug build, without MUMPS
            self.set_build_handler([
              'mkdir -p ${PREFIX}',
              './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=yes \
                --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
                --with-mpi-dir=${MPI_DIR} | tee out.txt',
              '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
              '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
            ])
          else:
            # release build without MUMPS
#              --with-blas-lapack-lib=${LAPACK_DIR}/lib/libopenblas.so\
            self.set_build_handler([
              'mkdir -p ${PREFIX}',
              './configure --prefix=${PREFIX} --with-shared-libraries=1 --with-debugging=no \
              --with-mpi-dir=${MPI_DIR}\
              COPTFLAGS=-O3\
              CXXOPTFLAGS=-O3\
              FOPTFLAGS=-O3 | tee out.txt',
            '$$(sed -n \'/Configure stage complete./{n;p;}\' out.txt) | tee out2.txt',
            '$$(sed -n \'/Now to install the libraries do:/{n;p;}\' out2.txt)',
            ])
          self.libs = ['petsc']
          
          self.number_output_lines = 3990
          
          res = super(PETSc, self).check(ctx, loc_callback=find_conf)
          self.check_required(res[0], ctx)
        
        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
