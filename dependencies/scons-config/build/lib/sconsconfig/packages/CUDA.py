import sys, os
from Package import Package

class CUDA(Package):

    def __init__(self, **kwargs):
        super(CUDA, self).__init__(**kwargs)
        self.libs=[['cuda', 'cudart']]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for CUDA ... ')
        self.check_options(env)

        res = super(CUDA, self).check(ctx)

        self.check_required(res[0])
        ctx.Result(res[0])
        return res[0]
