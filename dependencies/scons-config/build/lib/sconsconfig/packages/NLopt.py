import sys, os
from Package import Package


class NLopt(Package):

    def __init__(self, **kwargs):
        super(NLopt, self).__init__(**kwargs)
        self.libs=[['nlopt']]
        self.extra_libs=[]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <nlopt.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for NLopt ... ')
        self.check_options(env)

        res = super(NLopt, self).check(ctx)

        self.check_required(res[0])
        ctx.Result(res[0])
        return res[0]
