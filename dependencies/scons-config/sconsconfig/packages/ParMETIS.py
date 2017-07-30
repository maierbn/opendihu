import sys, os
from Package import Package

class ParMETIS(Package):

    def __init__(self, **kwargs):
        super(ParMETIS, self).__init__(**kwargs)
        self.libs=[['parmetis']]
        self.extra_libs=[['z', 'm', 'rt']]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <parmetis.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for ParMETIS ... ')
        self.check_options(env)

        res = super(ParMETIS, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
