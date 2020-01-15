import sys, os
from Package import Package

class boost(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'http://downloads.sourceforge.net/project/boost/boost/1.53.0/boost_1_53_0.tar.gz',
        }
        defaults.update(kwargs)
        super(boost, self).__init__(**defaults)
        self.ext = '.cc'
        self.sub_dirs = [
            ('include', 'lib'),
        ]
        self.headers = ['boost/optional.hpp']
        self.libs = [['boost_regex', 'boost_iostreams', 'boost_filesystem', 'boost_system']]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <boost/optional.hpp>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''
        self.set_build_handler([
            './bootstrap.sh',
            '!./b2 install --prefix=${PREFIX}'
        ])

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for boost ... ')
        self.check_options(env)

        res = super(boost, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
