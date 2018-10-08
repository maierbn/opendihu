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
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }
  
  LOG(TRACE) << "PythonFile::write ";

  // collect all available meshes
  std::set<std::string> meshNames;
  LoopOverTuple::loopCollectMeshNames<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), meshNames);
  
  // loop over meshes and create an output file for each
  for (std::string meshName : meshNames)
  {
    // setup name of file
    std::stringstream filenameStart;
    if (meshNames.size() == 1)
      filenameStart << this->filename_;
    else
     
      filenameStart << this->filename_ << "_" << meshName;
   
    // exelem file
    // determine file name
    std::stringstream s;
    s << filenameStart.str() << ".py";

    std::string filename = s.str();

    LOG(DEBUG) << "filename is [" << filename << "]";
    
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;
   
    // build python object for data
    PyObject *pyData = Python<typename DataType::FunctionSpace, typename DataType::OutputFieldVariables>::
      buildPyDataObject(data.getOutputFieldVariables(), meshName, timeStepNo, currentTime, this->onlyNodalValues_);
    //PyObject *pyData = PyDict_New();
    //PyDict_SetItemString(pyData, "a", PyLong_FromLong(5));
    //PyDict_SetItemString(pyData,"b", PyUnicode_FromString("hi"));

    if (VLOG_IS_ON(2))
    {
      VLOG(2) << "in python_file.tpp, data to write: " ;
      PythonUtility::printDict(pyData);
    }

    // open file, to see if directory needs to be created
    std::ofstream ofile = openFile(filename);
    if(ofile.is_open())
      ofile.close();

    // pickle is the python library to serialize objects
    bool usePickle = PythonUtility::getOptionBool(specificSettings_, "binary", false);

    std::string writeFlag = (usePickle? "wb" : "w");

    PyObject *file = openPythonFileStream(filename, writeFlag);

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
      outputPyObject(file, pyData);

      LOG(INFO) << (usePickle? "Binary" : "ASCII") << " file \"" << filename << "\" written.";
      LOG(DEBUG) << "PyFile_WriteObject done";
    }

    // decrement reference counters for python objects
    Py_XDECREF(pyData);
    Py_XDECREF(file);
    
  } // for meshName

#if PY_MAJOR_VERSION >= 3
  //PyGILState_Release(gilState);
#endif

  // for regular fixed also output stiffness matrix
  PythonStiffnessMatrixWriter<DataType>::writeNumpySolution(data, this->filename_);

}

};
