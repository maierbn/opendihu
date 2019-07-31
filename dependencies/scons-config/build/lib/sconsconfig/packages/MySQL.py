import sys, os
from Package import Package

##
##
##
class MySQL(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': '',
        }
        defaults.update(kwargs)
        super(MySQL, self).__init__(**defaults)
        self.ext = '.c'
        self.sub_dirs = [
            ('include/mysql', 'lib'),
            ('include/mysql', 'lib64'),
        ]
        self.headers = ['mysql.h']
        self.libs = ['mysqlclient']
        self.extra_libs = ['dl']
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <mysql.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for MySQL ... ')
        self.check_options(env)

        res = super(MySQL, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
