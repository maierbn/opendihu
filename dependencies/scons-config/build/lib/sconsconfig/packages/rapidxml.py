import sys, os
from Package import Package

def make_dir(x):
    if not os.path.exists(x):
        os.makedirs(x)

class rapidxml(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'http://downloads.sourceforge.net/project/rapidxml/rapidxml/rapidxml%201.13/rapidxml-1.13.zip'
        }
        defaults.update(kwargs)
        super(rapidxml, self).__init__(**defaults)
        self.ext = '.cc'
        self.sub_dirs = [('', ''), ('include', '')]
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <rapidxml.hpp>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''
        self.set_build_handler([
            (make_dir, '${PREFIX}'),
            'cp rapidxml.hpp rapidxml_iterators.hpp rapidxml_print.hpp rapidxml_utils.hpp ${PREFIX}',
        ])

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for rapidxml ... ')
        self.check_options(env)

        res = super(rapidxml, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
