import sys, os
from Package import Package


class GSL(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'http://gnu.mirror.uber.com.au/gsl/gsl-1.15.tar.gz',
        }
        defaults.update(kwargs)
        super(GSL, self).__init__(**defaults)
        self.libs=[
            ['gsl', 'gslcblas'],
            ['gsl'],
        ]
        self.extra_libs=[]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_version.h>
#include <gsl/gsl_integration.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

        # Setup the build handler. I'm going to assume this will work for all architectures.
        self.set_build_handler([
            './configure --prefix=${PREFIX}',
            'make',
            'make install'
        ])

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for GSL ... ')
        self.check_options(env)

        res = super(GSL, self).check(ctx)

        self.check_required(res[0])
        ctx.Result(res[0])
        return res[0]
