import sys, os
from Package import Package

##
##
##
class alsa(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': '',
        }
        defaults.update(kwargs)
        super(alsa, self).__init__(**defaults)
        self.ext = '.c'
        self.libs=[
            'asound',
        ]
        self.extra_libs=[
            [],
        ]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <alsa/asoundlib.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for alsa ... ')
        self.check_options(env)

        res = super(alsa, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
