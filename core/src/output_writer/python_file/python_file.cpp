#include "output_writer/python_file/python_file.h"

#include <Python.h>  // has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

PythonFile::PythonFile(DihuContext context, PythonConfig settings) : Generic(context, settings)
{
  onlyNodalValues_ = settings.getOptionBool("onlyNodalValues", true);
}

PyObject *PythonFile::openPythonFileStream(std::string filename, std::string writeFlag)
{
  //Note, this method already is called inside a critical section for the GIL
  
#if PY_MAJOR_VERSION >= 3
  //FILE *fileC = fopen(filename.c_str(), "w");
  //int fileDescriptorC = fileno(fileC);
  //PyObject *file = PyFile_FromFd(fileDescriptorC, NULL, "w", -1, NULL, NULL, NULL, true);

  //PyGILState_STATE gilState = PyGILState_Ensure();

  PyObject *ioModule = PyImport_ImportModule("io");

  PyObject *file = PyObject_CallMethod(ioModule, "open", "ss", filename.c_str(), writeFlag.c_str());
  Py_XDECREF(ioModule);

  LOG(DEBUG) << "file = io.open(" << filename << ", " << writeFlag << ")";
#else
  // python 2.7
  char filenameC[filename.size()+1];
  std::strcpy(filenameC, filename.c_str());

  char mode[writeFlag.size()+];
  std::strcpy(mode, writeFlag.c_str());

  PyObject *file = PyFile_FromString(filenameC, mode);
#endif

  if (!file)
  {
    LOG(ERROR) << "Could not open file \"" <<filename << "\" for output of python object.";
    return NULL;
  }

  return file;
}

void PythonFile::outputPyObject(PyObject *file, PyObject *pyData)
{

#if PY_MAJOR_VERSION >= 3
  // load json module
  static PyObject *jsonModule = NULL;
  if (jsonModule == NULL)
  {
    jsonModule = PyImport_ImportModule("json");
    LOG(DEBUG) << "import json";
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
    LOG(DEBUG) << "json.dump(data, file)";

    PyObject_CallMethod(file, "flush", NULL);
    PyObject_CallMethod(file, "close", NULL);
    LOG(DEBUG) << "file.close()";
  }
#else
  PyFile_WriteObject(pyData, file, 0);
#endif
}

}  // namespace
