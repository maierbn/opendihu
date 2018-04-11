import sys, os, multiprocessing
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
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          'cd ${SOURCE_DIR} && ./configure --prefix=${PREFIX} --enable-optimizations && make && make install',
          #'cd ${PREFIX} && make',
          #'mkdir -p ${PREFIX}/include && cp ${SOURCE_DIR}/*/*/*.h ${PREFIX}/include',
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
