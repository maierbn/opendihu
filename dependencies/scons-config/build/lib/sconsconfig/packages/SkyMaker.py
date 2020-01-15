import sys, os
from Package import Package

##
##
##
class SkyMaker(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': '',
        }
        defaults.update(kwargs)
        super(SkyMaker, self).__init__(**defaults)
        self.ext = '.c'
        self.libs = ['skymaker']
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <skymaker.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''
        # self.set_build_handler([
        #     './bootstrap.sh',
        #     '!./b2 install --prefix=${PREFIX}'
        # ])

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for SkyMaker ... ')
        self.check_options(env)

        res = super(SkyMaker, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
