import sys, os, multiprocessing, subprocess
from Package import Package

class Python(Package):

    def __init__(self, **kwargs):
        defaults = {
          #'download_url': 'https://www.python.org/ftp/python/2.7.14/Python-2.7.14.tgz',
          'download_url': 'https://www.python.org/ftp/python/2.7.13/Python-2.7.13.tgz',
          #https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz
        }
        defaults.update(kwargs)
        super(Python, self).__init__(**defaults)
        self.ext = '.cpp'
        #self.sub_dirs = [
        #    ('include/mysql', 'lib'),
        #    ('include/mysql', 'lib64'),
        #]
        #self.headers = ['mysql.h']
        #self.libs = ['mysqlclient']
        #self.extra_libs = ['lapack', 'blas']
        self.check_text = r'''
          #include <stdlib.h>
          #include <stdio.h>
          #include <Python.h>
          int main(int argc, char* argv[])
          {
            Py_Initialize();
            return EXIT_SUCCESS;
          }
        '''
    
        # Setup the build handler.
        self.libs = ["python2.7"]
        self.headers = ["Python.h"]
        
        # check configuration of gcc
        gcc_config = subprocess.check_output("gcc -v", shell=True)
        
        # if gcc was compiled such that -fuse-linker-plugin is available, compile with optimizations
        if "--enable-plugin" in gcc_config:        
          print "gcc has --enable-plugin, compile python with optimizations"
          self.set_build_handler([
            'mkdir -p ${PREFIX}',
            #'cd ${SOURCE_DIR} && ./configure --prefix=${PREFIX} --enable-optimizations && make && make install',
            'cd ${SOURCE_DIR} && ./configure --prefix=${PREFIX} && make && make install',
          ])
          self.number_output_lines = 7082
        else:       
          print "gcc has no --enable-plugin, compile python without optimizations"
          self.set_build_handler([
            'mkdir -p ${PREFIX}',
            'cd ${SOURCE_DIR} && ./configure --prefix=${PREFIX} && make && make install',
          ])
          self.number_output_lines = 7082

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Python ... ')
        self.check_options(env)

        res = super(Python, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
