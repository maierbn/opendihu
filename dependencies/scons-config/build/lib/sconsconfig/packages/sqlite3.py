import sys, os
from Package import Package

##
##
##
class sqlite3(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'http://github.com/furious-luke/sqlite3-ext/tarball/master',
        }
        defaults.update(kwargs)
        super(sqlite3, self).__init__(**defaults)
        self.ext = '.c'
        self.libs=[
            ['sqlite3'],
        ]
        self.extra_libs=[
            [],
        ]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <sqlite3.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

        # Setup the build handler. I'm going to assume this will work for all architectures.
        self.set_build_handler([
            './configure --prefix=${PREFIX}',
            'make install',
        ])

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for sqlite3 ... ')
        self.check_options(env)

        res = super(sqlite3, self).check(ctx)

        self.check_required(res[0])
        ctx.Result(res[0])
        return res[0]
