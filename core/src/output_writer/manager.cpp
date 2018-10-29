#include "output_writer/manager.h"

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "output_writer/python_callback/python_callback.h"
#include "output_writer/python_file/python_file.h"
#include "output_writer/paraview/paraview.h"
#include "output_writer/exfile/exfile.h"

namespace OutputWriter
{


void Manager::initialize(DihuContext context, PyObject *settings)
{
  outputWriter_.clear();

  VLOG(3) << "initializeOutputWriter, settings=" << settings;
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
      createOutputWriterFromSettings(context, writerSettings);
    }
  }
  else
  {
    LOG(DEBUG) << "Config does not contain \"OutputWriter\".";
  }
}

void Manager::createOutputWriterFromSettings(DihuContext context, PyObject *settings)
{

  if (PythonUtility::hasKey(settings, "format"))
  {
    // depending on type string create different OutputWriter object
    std::string typeString = PythonUtility::getOptionString(settings, "format", "none");
    if (typeString == "Paraview")
    {
      outputWriter_.push_back(std::make_shared<Paraview>(context, settings));
    }
    else if(typeString == "PythonCallback")
    {
      outputWriter_.push_back(std::make_shared<PythonCallback>(context, settings));
    }
    else if(typeString == "PythonFile")
    {
      outputWriter_.push_back(std::make_shared<PythonFile>(context, settings));
    }
    else if(typeString == "Exfile" || typeString == "ExFile")
    {
      outputWriter_.push_back(std::make_shared<Exfile>(context, settings));
    }
    else
    {
      LOG(WARNING) << "Unknown output writer type \"" << typeString<< "\". "
        << "Valid options are: \"Paraview\", \"PythonCallback\", \"PythonFile\", \"Exfile\"";
    }
  }
}

bool Manager::hasOutputWriters()
{
  return !outputWriter_.empty();
}

};
