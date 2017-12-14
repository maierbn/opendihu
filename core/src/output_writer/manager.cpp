#include "output_writer/manager.h"

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "output_writer/callback.h"
#include "output_writer/paraview.h"
#include "output_writer/python.h"

namespace OutputWriter
{
  
  
void Manager::initialize(PyObject *settings)
{
  outputWriter_.clear();
  
  LOG(DEBUG) << "initializeOutputWriter";
  PythonUtility::printDict(settings);
  
  if (PythonUtility::containsKey(settings, "OutputWriter"))
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
        outputWriter_.push_back(std::make_unique<Paraview>(settings));
      }
      else if(typeString == "Python")
      {
        outputWriter_.push_back(std::make_unique<Python>(settings));
      }
      else if(typeString == "Callback")
      {
        outputWriter_.push_back(std::make_unique<Callback>(settings));
      }
      else
      {
        LOG(WARNING) << "Unknown output writer type \""<<typeString<<"\".";
      }
    }
    else
    {
      LOG(WARNING) << "Output writer type is not a string";
    }
  }
  
  Py_CLEAR(key);
}

void Manager::writeOutput(Data::Data &problemData, int timeStepNo, double currentTime) const
{
  for(auto &outputWriter : outputWriter_)
  {
    outputWriter->write(problemData, timeStepNo, currentTime);
  }
}

};