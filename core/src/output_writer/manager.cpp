#include "output_writer/manager.h"

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "output_writer/python_callback/python_callback.h"
#include "output_writer/python_file/python_file.h"
#include "output_writer/paraview/paraview.h"
#include "output_writer/exfile/exfile.h"

namespace OutputWriter
{
  
  
void Manager::initialize(PyObject *settings)
{
  outputWriter_.clear();
  
  LOG(DEBUG) << "initializeOutputWriter, settings=" << settings;
  //PythonUtility::printDict(settings);
  
  if (PythonUtility::hasKey(settings, "OutputWriter"))
  {
    // get the first value from the list
    PyObject *writerSettings = PythonUtility::getOptionListBegin<PyObject *>(settings, "OutputWriter");
  
    // loop over other values
    for (;
        !PythonUtility::getOptionListEnd(settings, "OutputWriter");
        PythonUtility::getOptionListNext<PyObject *>(settings, "OutputWriter", writerSettings))
    {
      LOG(DEBUG) << "parse outputWriter";
      createOutputWriterFromSettings(writerSettings);
    }
  }
  else
  {
    LOG(DEBUG) << "config does not contain \"OutputWriter\".";
  }
}

void Manager::createOutputWriterFromSettings(PyObject *settings)
{
  PyObject *key = PyString_FromString("format");
  if (PyDict_Contains(settings, key))
  {
    PyObject *type = PyDict_GetItem(settings, key);
    if (PyString_Check(type))
    {
      // depending on type string create different OutputWriter object
      std::string typeString = PyString_AsString(type);
      if (typeString == "Paraview")
      {
        outputWriter_.push_back(std::make_shared<Paraview>(settings));
      }
      else if(typeString == "PythonCallback")
      {
        outputWriter_.push_back(std::make_shared<PythonCallback>(settings));
      }
      else if(typeString == "PythonFile")
      {
        outputWriter_.push_back(std::make_shared<PythonFile>(settings));
      }
      else if(typeString == "Exfile")
      {
        outputWriter_.push_back(std::make_shared<Exfile>(settings));
      }
      else
      {
        LOG(WARNING) << "Unknown output writer type \""<<typeString<<"\". " 
          << "Valid options are: \"Paraview\", \"PythonCallback\", \"PythonFile\", \"Exfile\"";
      }
    }
    else
    {
      LOG(WARNING) << "Output writer type is not a string";
    }
  }
  
  Py_CLEAR(key);
}


};