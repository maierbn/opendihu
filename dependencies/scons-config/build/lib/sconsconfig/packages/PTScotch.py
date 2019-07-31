import sys, os
from Package import Package

class PTScotch(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://gforge.inria.fr/frs/download.php/28977/scotch_5.1.12b.tar.gz'
        }
        defaults.update(kwargs)
        super(PTScotch, self).__init__(**defaults)
        self.libs=[['ptscotch', 'ptscotcherr', 'ptscotcherrexit', 'scotch', 'scotcherr', 'scotcherrexit']]
        self.extra_libs=[['z', 'm'],
                         ['z', 'm', 'rt']]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <ptscotch.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

        # Setup the build handlers. Thanks to the way PTScotch is setup we will need a unique one of these
        # for every bloody architecture under the sun.
        replace_list = [
            'sed -i \'s/-DPTHREAD_COMMON//g\' Makefile.inc', # remove PTHREAD
            'sed -i \'s/-DSCOTCH_PTHREAD//g\' Makefile.inc', # remove PTHREAD
            'sed -i \'s!mpicc!${MPI_DIR}/bin/mpicc!g\' Makefile.inc', # replace mpicc with location
            'sed -i \'s!gcc!${MPI_DIR}/bin/mpicc!g\' Makefile.inc', # replace gcc with mpicc
            'sed -i \'s!/usr/local!${PREFIX}!g\' Makefile', # replace with prefix
        ]
        make_list = [
            'make',
            'make ptscotch',
            (lambda x: os.path.exists(x) or os.mkdir(x), '${PREFIX}'),
            'make install',
            (os.chdir, '..'), # must back out!
        ]
        self.set_build_handler([
            (os.chdir, 'src'),
            'cp Make.inc/Makefile.inc.x86-64_pc_linux2 Makefile.inc',
        ] + replace_list + make_list)
        self.set_build_handler([
            (os.chdir, 'src'),
            'cp Make.inc/Makefile.inc.i686_mac_darwin8 Makefile.inc',
        ] + replace_list + make_list, 'Darwin')

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for PTScotch ... ')
        self.check_options(env)

        res = super(PTScotch, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
