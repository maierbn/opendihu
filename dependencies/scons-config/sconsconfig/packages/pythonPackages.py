import sys, os
from .Package import Package

check_text = r'''

#include <Python.h>  // this has to be the first included header

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

int main(int argc, char* argv[])
{
  std::string pythonSearchPath = "${DEPENDENCIES_DIR}/python/install";
  
  std::cout << "pythonSearchPath: [" << pythonSearchPath << "]" << std::endl;
  
  //std::string pythonSearchPath = std::string("/store/software/opendihu/dependencies/python/install");
  const wchar_t *pythonSearchPathWChar = Py_DecodeLocale(pythonSearchPath.c_str(), NULL);
  Py_SetPythonHome((wchar_t *)pythonSearchPathWChar);

  Py_Initialize();
  
  PyEval_InitThreads();
  Py_SetStandardStreamEncoding(NULL, NULL);

  // check if numpy module could be loaded
  PyObject *numpyModule = PyImport_ImportModule("numpy");
  if (numpyModule == NULL)
  {
    std::cout << "Failed to import numpy." << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
'''

class pythonPackages(Package):
  
    def __init__(self, **kwargs):
        defaults = {
          'download_url': 'https://files.pythonhosted.org/packages/50/2f/618c5b6804e6dda90f024f0bc2d14ffc3db00221a818ee35da478427015d/svg.path-3.0-py2.py3-none-any.whl',
        }
        defaults.update(kwargs)
        super(pythonPackages, self).__init__(**defaults)
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
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install --upgrade pip',
              '$${DEPENDENCIES_DIR}/python/install/bin/python3 -m pip install ${DEPENDENCIES_DIR}/pythonpackages/*.whl --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
        else :
          # Setup the build handler.
          self.set_build_handler([
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install --upgrade pip',
              '$${DEPENDENCIES_DIR}/python/install/bin/python3 -m pip install numpy matplotlib scipy numpy-stl svg.path triangle geomdl pymp vtk --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
        
        self.number_output_lines = 13780
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Python packages...')
        self.check_options(env)

        res = super(pythonPackages, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
