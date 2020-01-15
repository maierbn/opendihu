import sys, os
from Package import Package

##
##
##
class libsamplerate(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': '',
        }
        defaults.update(kwargs)
        super(libsamplerate, self).__init__(**defaults)
        self.ext = '.c'
        self.libs=[
            'samplerate',
        ]
        self.extra_libs=[
            [],
        ]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <samplerate.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for libsamplerate ... ')
        self.check_options(env)

        res = super(libsamplerate, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
