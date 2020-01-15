import sys, os
from Package import Package

##
##
##
class clock(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': '',
        }
        defaults.update(kwargs)
        super(clock, self).__init__(**defaults)
        self.ext = '.c'
        self.libs=[
            [],
            ['rt'],
        ]
        self.extra_libs=[]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
int main(int argc, char* argv[]) {
   struct timespec ts;
   clock_gettime( CLOCK_MONOTONIC, &ts );
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for clock ... ')
        self.check_options(env)

        res = super(clock, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
