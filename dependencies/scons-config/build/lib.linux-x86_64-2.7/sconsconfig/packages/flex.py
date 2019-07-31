import sys, os, multiprocessing
from Package import Package
import subprocess

class flex(Package):

  def __init__(self, **kwargs):
    defaults = {
      'download_url': 'https://github.com/westes/flex/releases/download/v2.6.1/flex-2.6.1.tar.gz',
    }
    defaults.update(kwargs)
    super(flex, self).__init__(**defaults)
    self.ext = '.cpp'
    #self.sub_dirs = [
    #    ('include/mysql', 'lib'),
    #    ('include/mysql', 'lib64'),
    #]
    #self.headers = ['mysql.h']
    #self.libs = ['mysqlclient']
    #self.extra_libs = ['flex', 'blas']
    #self.set_rpath = False    # do not use dynamic linkage
    self.check_text = r'''
#include <cstdlib>
#include <iostream>
#include <FlexLexer.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''
    self.set_build_handler([
      'PATH=${PATH}:${DEPENDENCIES_DIR}/bison/install/bin ./autogen.sh',
      'PATH=${PATH}:${DEPENDENCIES_DIR}/bison/install/bin ./configure --prefix=${PREFIX} ',
      'PATH=${PATH}:${DEPENDENCIES_DIR}/bison/install/bin make all',
      'make install'
    ])
    self.number_output_lines = 577

    self.libs = ["fl", "fl_pic"]
    self.headers = ["FlexLexer.h"]

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for flex ... ')
    self.check_options(env)

    res = super(flex, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
