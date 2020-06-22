import sys, os, multiprocessing
from .Package import Package

class libxml2(Package):

    def __init__(self, **kwargs):
        defaults = {
            #'download_url': 'https://git.gnome.org/browse/libxml2/snapshot/libxml2-2.9.7.zip',
            'download_url': 'ftp://xmlsoft.org/libxml2/libxml2-git-snapshot.tar.gz',
        }
        defaults.update(kwargs)
        super(libxml2, self).__init__(**defaults)
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <libxml/xmlmodule.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

        # get number of available processors
        p = multiprocessing.cpu_count()

        # Setup the build handler.
        self.set_build_handler([
          "cd ${SOURCE_DIR} && ./configure --prefix=${PREFIX} --without-python && make && make install",
          "mv ${PREFIX}/include/libxml2/* ${PREFIX}/include && rmdir ${PREFIX}/include/libxml2"
        ])
        self.number_output_lines = 521
          
        self.libs = ["xml2"]
        self.headers = ["libxml/xmlmodule.h"]

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for libxml2 ...       ')
        self.check_options(env)

        res = super(libxml2, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
