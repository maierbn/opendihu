import sys, os
from Package import Package

class AntTweakBar(Package):

    def __init__(self, **kwargs):
        defaults = {
        }
        defaults.update(kwargs)
        super(AntTweakBar, self).__init__(**defaults)
        self.libs=[
            ['AntTweakBar'],
        ]
        self.extra_libs=[]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <AntTweakBar.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for AntTweakBar ... ')
        self.check_options(env)

        res = super(AntTweakBar, self).check(ctx)

        self.check_required(res[0])
        ctx.Result(res[0])
        return res[0]
