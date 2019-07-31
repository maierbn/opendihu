import sys, os
from Package import Package


class FFTW(Package):

    def __init__(self, **kwargs):
        super(FFTW, self).__init__(**kwargs)
        self.libs=[['srfftw_mpi', 'sfftw_mpi', 'srfftw', 'sfftw']]
        self.extra_libs=[]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <srfftw_mpi.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for FFTW ... ')
        self.check_options(env)

        res = super(FFTW, self).check(ctx)

        self.check_required(res[0])
        ctx.Result(res[0])
        return res[0]
