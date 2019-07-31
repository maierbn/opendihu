import sys, os
from Package import Package

##
##
##
class PostgreSQL(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': '',
        }
        defaults.update(kwargs)
        super(PostgreSQL, self).__init__(**defaults)
        self.ext = '.c'
        # self.sub_dirs = [
        #     ('include/mysql', 'lib'),
        #     ('include/mysql', 'lib64'),
        # ]
        self.headers = ['libpq-fe.h']
        self.libs = ['pq']
        # self.extra_libs = ['dl']
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <libpq-fe.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for PostgreSQL ... ')
        self.check_options(env)

        res = super(PostgreSQL, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
