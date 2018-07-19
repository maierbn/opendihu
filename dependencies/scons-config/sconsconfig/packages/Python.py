import sys, os, multiprocessing, subprocess
from Package import Package

#
# Python requires bzip2 (for matplotlib to work)
#
class Python(Package):

    def __init__(self, **kwargs):
        defaults = {
          #'download_url': 'https://www.python.org/ftp/python/2.7.14/Python-2.7.14.tgz',
          #'download_url': 'https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz',
          #'download_url': 'https://www.python.org/ftp/python/3.4.6/Python-3.4.6.tgz'
          'download_url': 'https://www.python.org/ftp/python/3.6.5/Python-3.6.5.tgz'
          #'download_url': 'https://www.python.org/ftp/python/3.5.3/Python-3.5.3.tgz'
            
        }
        defaults.update(kwargs)
        super(Python, self).__init__(**defaults)
        self.ext = '.cpp'
        self.sub_dirs = [
            ('include', 'lib'),
        #    ('include/mysql', 'lib64'),
        ]
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
        self.libs = ["python3.6m"]
        self.headers = ["Python.h"]
        
        # check configuration of gcc
        gcc_config = subprocess.check_output(["gcc", "-v"], stderr=subprocess.STDOUT)
        
        # extract and output gcc version
        pos1 = gcc_config.find("gcc version")
        pos2 = pos1+11+gcc_config[pos1+11:].find(" ")
        pos3 = pos2+gcc_config[pos2+1:].find(" ")
        #print("GCC version: {}".format(gcc_config[pos2+1:pos3+1]))
        #gcc_config = ""
        
        # if gcc was compiled such that -fuse-linker-plugin is available, compile with optimizations
        if "--enable-plugin" in gcc_config:        
          #print("gcc has --enable-plugin, compile python with optimizations")
          self.set_build_handler([
            'mkdir -p ${PREFIX}',
            'cd ${SOURCE_DIR} && ./configure --enable-shared --enable-optimizations --prefix=${PREFIX} \
              LDFLAGS="-Wl,--rpath=${PREFIX}/lib -L${DEPENDENCIES_DIR}/bzip2/install/lib" \
              CFLAGS="-I${DEPENDENCIES_DIR}/bzip2/install/include" \
              && make && make install',
            '$export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PREFIX}/lib',
            'cd ${PREFIX}/include && echo "#define PYTHON_HOME_DIRECTORY \\"${PREFIX}\\"\n" > python_home.h',
          ])
          self.number_output_lines = 9823
        else:       
          print("gcc has no --enable-plugin, compile python without optimizations")
          self.set_build_handler([
            'mkdir -p ${PREFIX}',
            'cd ${SOURCE_DIR} && ./configure --enable-shared --prefix=${PREFIX} \
              LDFLAGS="-Wl,--rpath=${PREFIX}/lib -L${DEPENDENCIES_DIR}/bzip2/install/lib" \
              CFLAGS="-I${DEPENDENCIES_DIR}/bzip2/install/include" \
              && make && make install',
            '$export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PREFIX}/lib',
            'cd ${PREFIX}/include && echo "#define PYTHON_HOME_DIRECTORY \\"${PREFIX}\\"\n" > python_home.h',
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
