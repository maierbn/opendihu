import sys, os
from Package import Package

class libusb(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': '',
        }
        defaults.update(kwargs)
        super(libusb, self).__init__(**defaults)
        self.ext = '.c'
        self.libs = ['usb-1.0']
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <libusb-1.0/libusb.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for libusb ... ')
        self.check_options(env)

        res = super(libusb, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
