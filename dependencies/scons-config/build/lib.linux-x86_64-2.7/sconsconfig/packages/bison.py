import sys, os, multiprocessing
from Package import Package
import subprocess

class bison(Package):

  def __init__(self, **kwargs):
    defaults = {
      'download_url': 'http://ftp.gnu.org/gnu/bison/bison-3.1.tar.gz',
    }
    defaults.update(kwargs)
    super(bison, self).__init__(**defaults)
    #self.ext = '.c'
    #self.sub_dirs = [
    #    ('include/mysql', 'lib'),
    #    ('include/mysql', 'lib64'),
    #]
    #self.headers = ['mysql.h']
    #self.libs = ['mysqlclient']
    #self.extra_libs = ['bison', 'blas']
    #self.set_rpath = False    # do not use dynamic linkage
    self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''
    self.set_build_handler([
      './configure --prefix=${PREFIX} ',
      'make all',
      'make install',
      'export PATH=$PATH:${PREFIX}/bin'
    ])
    self.number_output_lines = 743

    #self.libs = ["openblas"]
    #self.headers = ["bisone.h"]

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for bison ... ')
    self.check_options(env)

    res = super(bison, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
