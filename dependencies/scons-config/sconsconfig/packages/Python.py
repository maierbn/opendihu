import sys, os, multiprocessing, subprocess, socket
from .Package import Package

#
# Python requires bzip2 (for matplotlib to work)
#
class Python(Package):

  def __init__(self, **kwargs):
    defaults = {
      #'download_url': 'https://www.python.org/ftp/python/2.7.14/Python-2.7.14.tgz',
      #'download_url': 'https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz',
      #'download_url': 'https://www.python.org/ftp/python/3.4.6/Python-3.4.6.tgz'
      #'download_url': 'https://www.python.org/ftp/python/3.5.3/Python-3.5.3.tgz'
      #'download_url': 'https://www.python.org/ftp/python/3.6.5/Python-3.6.5.tgz'
      #'download_url': 'https://github.com/stackless-dev/stackless/archive/3.6-slp.zip'  # stackless-python
      'download_url': 'https://www.python.org/ftp/python/3.9.0/Python-3.9.0.tgz'
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

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for Python 3.9 ...   ')
    
    # python 3.9
    self.libs = ["python3.9"]
    self.headers = ["Python.h"]
 
    self.set_build_handler([
      'mkdir -p ${PREFIX}',
      'cd ${SOURCE_DIR} && chmod +x ./configure && ./configure --enable-shared --prefix=${PREFIX} \
        LDFLAGS="-Wl,--rpath=${PREFIX}/lib -L${DEPENDENCIES_DIR}/bzip2/install/lib -L${DEPENDENCIES_DIR}/zlib/install/lib" \
        CPPFLAGS="-I${DEPENDENCIES_DIR}/bzip2/install/include -I${DEPENDENCIES_DIR}/zlib/install/include" CXX=g++ CC=gcc \
        && make && make install',
      '$export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PREFIX}/lib',
      'cd ${PREFIX}/include && echo "#define PYTHON_HOME_DIRECTORY \\"${PREFIX}\\"\n" > python_home.h',
    ])
    self.number_output_lines = 11540

    self.check_options(env)
    res = super(Python, self).check(ctx)

    # if installation failed with the current command, retry with different options 
    if not res[0]:
      ctx.Log('Retry with python 3.6\n')
      ctx.Message('Retry with python 3.6 ...')

      # python 3.6
      self.libs = ["python3.6m"]
      self.download_url = 'https://www.python.org/ftp/python/3.6.5/Python-3.6.5.tgz'

      self.check_options(env)
      res = super(Python, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
