#include "output_writer/python_file/python_file.h"

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <cstdio>

#include "easylogging++.h"
#include "output_writer/python/python.h"
#include "output_writer/python_file/python_stiffness_matrix_writer.h"

namespace OutputWriter
{

template<typename DataType>
void PythonFile::write(DataType& data, int timeStepNo, double currentTime)
{ 
  LOG(TRACE) << "PythonFile::write ";
 
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }
  
  // build python object for data
  PyObject *pyData = Python<typename DataType::BasisOnMesh, typename DataType::OutputFieldVariables>::
    buildPyDataObject(data.getOutputFieldVariables(), timeStepNo, currentTime, this->onlyNodalValues_);
  //PyObject *pyData = PyDict_New();
  //PyDict_SetItemString(pyData, "a", PyLong_FromLong(5));
  //PyDict_SetItemString(pyData,"b", PyUnicode_FromString("hi"));
  
  //PythonUtility::printDict(pyData);
  
  // determine file name
  std::stringstream s;
  s<<this->filename_<<".py";
  std::string filename = s.str();
  
  LOG(DEBUG) << "filename is [" << filename << "]";
  
  // open file, to see if directory needs to be created
  std::ofstream ofile = openFile(filename);
  if(ofile.is_open())
    ofile.close();
  
#if PY_MAJOR_VERSION >= 3  
  //FILE *fileC = fopen(filename.c_str(), "w");
  //int fileDescriptorC = fileno(fileC);
  //PyObject *file = PyFile_FromFd(fileDescriptorC, NULL, "w", -1, NULL, NULL, NULL, true);
  
  //PyGILState_STATE gilState = PyGILState_Ensure();

  PyObject *ioModule = PyImport_ImportModule("io");

  PyObject *file = PyObject_CallMethod(ioModule, "open", "ss", filename.c_str(), "w");
  Py_DECREF(ioModule);
#else 
  // python 2.7
  char filenameC[filename.size()+1];
  std::strcpy(filenameC, filename.c_str());
  
  char mode[3];
  std::strcpy(mode, "wb");

  PyObject *file = PyFile_FromString(filenameC, mode);
#endif
 
  if (!file) 
  {
    LOG(ERROR) << "Could not open file \""<<filename<<"\" for output of python object.";
    return;
  }
  
  // pickle is the python library to serialize objects
  bool usePickle = PythonUtility::getOptionBool(specificSettings_, "binary", false);
  
  if (usePickle)
  {
    // load pickle module if is was not loaded in an earlier call
    static PyObject *pickleModule = NULL;
    if (pickleModule == NULL)
    {
      pickleModule = PyImport_ImportModuleNoBlock("pickle");
    }
    if (pickleModule == NULL)
    {
      LOG(ERROR) << "Could not import pickle module";
    }
    else
    {
     /*
      std::string methodNameStr("dump");
      char methodName[methodNameStr.size()+1];
      std::strcpy(methodName, methodNameStr.c_str());
      
      std::string formatStr("(O O i)");
      char format[formatStr.size()+1];
      std::strcpy(format, formatStr.c_str());
      */
      PyObject *pickle = PyObject_CallMethod(pickleModule, "dump", "(O O i)", pyData, file, 1);
      Py_XDECREF(pickle);
    }
  }
  else
  {
    // output python object
   
#if PY_MAJOR_VERSION >= 3 
    // load json module
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
#else
    PyFile_WriteObject(pyData, file, 0);
#endif
    
    LOG(INFO) << (usePickle? "Binary" : "ASCII") << " file \"" << filename << "\" written.";
    LOG(DEBUG) << "PyFile_WriteObject done";
  }
  
  // decrement reference counters for python objects
  Py_XDECREF(pyData);
  Py_XDECREF(file);
  
#if PY_MAJOR_VERSION >= 3  
  //PyGILState_Release(gilState);
#endif
  
  LOG(DEBUG) << "writeNumpySolution";
  
  // for regular fixed also output stiffness matrix 
  PythonStiffnessMatrixWriter<DataType>::writeNumpySolution(data, this->filename_);
}

};