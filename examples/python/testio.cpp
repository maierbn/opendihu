#include "opendihu.h"

#include <iostream>

int main(int argc, char* argv[])
{
  Py_Initialize();
  
  // test output
  const char *version = Py_GetVersion();
  std::cout << "python version: " << version << std::endl;
    
  PyObject *pyData = PyDict_New();
  PyDict_SetItemString(pyData, "a", PyLong_FromLong(5));
  PyDict_SetItemString(pyData, "b", PyUnicode_FromString("hi"));
  
  PyObject *ioModule = PyImport_ImportModule("io");

  std::string filename = "out_testio.py";
  PyObject *file = PyObject_CallMethod(ioModule, "open", "ss", filename.c_str(), "w");
  Py_DECREF(ioModule);
  
  static PyObject *jsonModule = NULL;
  if (jsonModule == NULL)
  {
    jsonModule = PyImport_ImportModule("json");
  }
  if (jsonModule == NULL)
  {
    LOG(ERROR) << "Could not import json module";
  }
  else
  {
    // convert data object to string representation and save in file
    //PyObject *pyDataJson = 
    PyObject_CallMethod(jsonModule, "dump", "O, O", pyData, file);
    
    PyObject_CallMethod(file, "flush", NULL);
    PyObject_CallMethod(file, "close", NULL);
  }
  
  return EXIT_SUCCESS;
}
