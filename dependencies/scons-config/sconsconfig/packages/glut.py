import sys, os
from Package import Package

class glut(Package):

    def __init__(self, **kwargs):
        defaults = {
        }
        defaults.update(kwargs)
        super(glut, self).__init__(**defaults)
        self.libs=[
            ['GL', 'GLU', 'glut'],
        ]
        self.extra_libs=[]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <glut.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for glut ... ')
        self.check_options(env)

        res = super(glut, self).check(ctx)

        self.check_required(res[0])
        ctx.Result(res[0])
        return res[0]
