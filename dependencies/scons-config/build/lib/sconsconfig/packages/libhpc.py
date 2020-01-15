import sys, os
from Package import Package
import sconsconfig as config

##
##
##
class libhpc(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'http://github.com/furious-luke/libhpc/tarball/master',
        }
        defaults.update(kwargs)
        super(libhpc, self).__init__(**defaults)
        self.ext = '.cc'
        self.libs=[
            'hpc',
            {'libraries': 'hpc', 'prepend': False},
        ]
        self.extra_libs=[
            [],
        ]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <libhpc/libhpc.hh>
int main(int argc, char* argv[]) {
   int rank;
   MPI_Init(&argc, &argv);
   printf("%d\n", MPI_VERSION);
   printf("%d\n", MPI_SUBVERSION);

   /* This will catch OpenMPI libraries. */
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   MPI_Finalize();
   return EXIT_SUCCESS;
}
'''

        # MPI on a cluster usually won't run properly, so don't
        # try to.
        self.run = False

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for libhpc ... ')
        self.check_options(env)

        # Check for a bunch of things we need for auto building.
        cmd = 'scons PREFIX=${PREFIX}'
        pkg = config.package(config.packages.boost)
        if pkg and pkg.found and pkg.base_dir:
            cmd += ' BOOST_DIR=' + pkg.base_dir
        pkg = config.package(config.packages.MPI)
        if pkg and pkg.found and pkg.base_dir:
            cmd += ' MPI_DIR=' + pkg.base_dir
        pkg = config.package(config.packages.HDF5)
        if pkg and pkg.found and pkg.base_dir:
            cmd += ' HDF5_DIR=' + pkg.base_dir
        pkg = config.package(config.packages.pugixml)
        if pkg and pkg.found and pkg.base_dir:
            cmd += ' PUGIXML_DIR=' + pkg.base_dir
        if env['WITH_OPENMP']:
            cmd += ' WITH_OPENMP=yes'
        cmd += ' install'

        # Setup the build handler. I'm going to assume this will work for all architectures.
        self.set_build_handler([cmd])

        res = super(libhpc, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
