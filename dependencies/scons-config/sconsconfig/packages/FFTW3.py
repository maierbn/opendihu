import sys, os
from Package import Package


class FFTW3(Package):

    def __init__(self, **kwargs):
        use_float = kwargs.pop('use_float', False)
        super(FFTW3, self).__init__(**kwargs)
        if use_float:
            self.libs = ['fftw3f']
        else:
            self.libs = ['fftw3']
        self.extra_libs = []
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for FFTW3 ... ')
        self.check_options(env)

        res = super(FFTW3, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
