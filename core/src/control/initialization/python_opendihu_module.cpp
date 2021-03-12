#include "control/initialization/python_opendihu_module.h"

#include <functional>
#include <iostream>
#include <string>
#include <Python.h>

namespace PythonOpendihu
{

PyModuleDef opendihuModule =
{
    PyModuleDef_HEAD_INIT,
    "opendihu", 0, -1, 0,
};

// Internal state
std::string programName;      //< the name of the executable, as given in argv[0]
std::string versionText;
std::string metaText;

PyMODINIT_FUNC PyInit_opendihu(void)
{
  // create the python object of the new "opendihu" module
  PyObject* module = PyModule_Create(&opendihuModule);

  // if the module could be created successfully, add string constants
  if (module)
  {
    PyModule_AddStringConstant(module, "program_name", programName.c_str());
    PyModule_AddStringConstant(module, "version", versionText.c_str());
    PyModule_AddStringConstant(module, "meta", metaText.c_str());
  }
  return module;
}

} // namespace opendihu

