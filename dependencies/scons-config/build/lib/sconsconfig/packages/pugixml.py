import sys, os
from Package import Package

def make_dir(x):
    if not os.path.exists(x):
        os.makedirs(x)

class pugixml(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'http://pugixml.googlecode.com/files/pugixml-1.2.tar.gz'
        }
        defaults.update(kwargs)
        super(pugixml, self).__init__(**defaults)
        self.ext = '.cc'
        # self.sub_dirs = [('', ''), ('include', '')]
        self.libs = ['pugixml']
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <pugixml.hpp>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''
        self.set_build_handler([
            (make_dir, '${PREFIX}'),
            (make_dir, '${PREFIX}/include'),
            (make_dir, '${PREFIX}/lib'),
            'cmake scripts -DBUILD_SHARED_LIBS=yes',
            'make',
            'cp src/pugixml.hpp src/pugiconfig.hpp ${PREFIX}/include',
            'cp libpugixml.so.1.2 ${PREFIX}/lib',
            'ln -s ${PREFIX}/lib/libpugixml.so.1.2 ${PREFIX}/lib/libpugixml.so',
        ])

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for pugixml ... ')
        self.check_options(env)

        self.need_cmake(env)
        res = super(pugixml, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
