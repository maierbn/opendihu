import sys, os
from distutils import sysconfig
from .Package import Package

check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <Python.h>
#include <iostream>

int main(int argc, char* argv[]) {

  Py_Initialize();
  
  std::string text = R"(
import svg.path
)";
  int ret = 0;
  try
  {
    std::cout << "run test program" << std::endl;
    ret = PyRun_SimpleString(text.c_str());
  }
  catch(...)
  {
    std::cout << "test program failed, ret=" << ret << std::endl;
  }
  std::cout << std::string(80, '-') << std::endl;
  
  // if there was an error in the python code
  if (ret != 0)
  {
    if (PyErr_Occurred())
    {
      // print error message and exit
      PyErr_Print();
      std::cout << "An error occurred in the python test script." << std::endl;
      exit(EXIT_FAILURE);
    }
    exit(EXIT_FAILURE);
  }
   return EXIT_SUCCESS;
}
'''

class SvgPath(Package):
  
    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://files.pythonhosted.org/packages/50/2f/618c5b6804e6dda90f024f0bc2d14ffc3db00221a818ee35da478427015d/svg.path-3.0-py2.py3-none-any.whl'
        }
        defaults.update(kwargs)
        super(SvgPath, self).__init__(**defaults)
        self.ext = '.cpp'
        self.sub_dirs = [
            ('include', 'lib'),
        ]
        self.headers = []
        self.libs = []
        self.check_text = check_text
        self.static = False
        
        if os.environ.get("PE_ENV") is not None:
          # Setup the build handler.
          self.set_build_handler([
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install ${PREFIX}/../svg.path-3.0-py2.py3-none-any.whl --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
        else :
          # Setup the build handler.
          self.set_build_handler([
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install svg.path --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
        
        # SvgPath is installed in the directory tree of python, under lib/python3.6/site-packages. It is done using pip.

        self.number_output_lines = 13780
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for svg.path ...      ')
        self.check_options(env)

        res = super(SvgPath, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
