#include "opendihu.h"

#include <iostream>

int main(int argc, char* argv[])
{
  int ret = 0;
  
  //wchar_t *programNameWChar = Py_DecodeLocale(argv[0], NULL);
  //Py_SetProgramName(programNameWChar);  /* optional but recommended */
    
  //std::string pythonSearchPath = std::string("/store/software/opendihu/dependencies/python/install");
  //const wchar_t *pythonSearchPathWChar = Py_DecodeLocale(pythonSearchPath.c_str(), NULL);
  //Py_SetPythonHome((wchar_t *)pythonSearchPathWChar);
  
  //std::string pythonPath("/store/software/opendihu/dependencies/python/install/lib/python3.6");
  //const wchar_t *pythonPathWChar = Py_DecodeLocale(pythonPath.c_str(), NULL);
  //Py_SetPath((wchar_t *)pythonPathWChar);
  
  //dlopen("/store/software/opendihu/dependencies/python/install/lib/libpython3.6m.so", RTLD_LAZY | RTLD_GLOBAL);
  
  Py_Initialize();
  
  wchar_t *pathWChar = Py_GetPath();
  char *path = Py_EncodeLocale(pathWChar, NULL);
  std::cout << "python path: " << path << std::endl;
  
  // test output
  PyObject *ioModule = PyImport_ImportModule("io");
  if (ioModule == NULL)
    std::cout << "could not load io module" << std::endl;
/*
/store/software/opencmiss/own_scripts
/store/software/opendihu/scripts
/store/software/opendihu/dependencies/python/install/lib/python36.zip
/store/software/opendihu/dependencies/python/install/lib/python3.6
/store/software/opendihu/dependencies/python/install/lib/python3.6/lib-dynload
/store/software/opendihu/dependencies/python/install/lib/python3.6/site-packages
/store/software/opendihu/dependencies/python/install/lib/python3.6/site-packages/Cython-0.27.2-py3.6-linux-x86_64.egg
/store/software/opendihu/dependencies/python/install/lib/python3.6/site-packages/numpy-1.15.0.dev0+unknown-py3.6-linux-x86_64.egg

/store/software/opencmiss/own_scripts
/store/software/opendihu/scripts
/store/software/opendihu/dependencies/python/install/lib/python3.6/site-packages
/store/software/opendihu/dependencies/python/install/lib/python36.zip
/store/software/opendihu/dependencies/python/install/lib/python3.6
/store/software/opendihu/dependencies/python/install/lib/python3.6/lib-dynload
/store/software/opendihu/dependencies/python/install/lib/python3.6/site-packages/Cython-0.27.2-py3.6-linux-x86_64.egg
/store/software/opendihu/dependencies/python/install/lib/python3.6/site-packages/numpy-1.15.0.dev0+unknown-py3.6-linux-x86_64.egg
*/
  
  //Py_SetStandardStreamEncoding(NULL, NULL);
  
  //wchar_t **argvWChar = new wchar_t *[argc];
  //for (int i=0; i<argc; i++)
  //{
  //  argvWChar[i] = Py_DecodeLocale(argv[i], NULL);
  //}

  //PySys_SetArgvEx(argc, argvWChar, 1);
  //PyEval_InitThreads();

  const char *version = Py_GetVersion();
  std::cout << "python version: " << version << std::endl;
    
  std::cout << std::string(80, '-') << std::endl;
  
  std::string text = R"(
print("hello")

import sys
print("pythonpath:")
for p in sys.path:
  print(p)

import numpy as np
print(np.pi)

import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
print("installed packages: ",installed_packages_list)

)";
  //int ret = 0;
  try
  {
    ret = PyRun_SimpleString(text.c_str());
  }
  catch(...)
  {
  }
  std::cout << std::string(80, '-') << std::endl;
  
  // if there was an error in the python code
  if (ret != 0)
  {
    if (PyErr_Occurred())
    {
      // print error message and exit
      PyErr_Print();
      exit(EXIT_FAILURE);
    }
    exit(EXIT_FAILURE);
  }
  
  return EXIT_SUCCESS;
}
