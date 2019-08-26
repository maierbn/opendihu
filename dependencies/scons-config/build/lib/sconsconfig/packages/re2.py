import sys, os
from Package import Package

class re2(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': '',
        }
        defaults.update(kwargs)
        super(re2, self).__init__(**defaults)
        self.ext = '.cc'
        self.libs = ['re2']
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <re2/re2.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for re2 ... ')
        self.check_options(env)

        res = super(re2, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
