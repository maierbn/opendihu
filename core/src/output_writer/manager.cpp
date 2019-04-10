#include "output_writer/manager.h"

#include "easylogging++.h"

#include "control/python_home.h"
#include "utility/python_utility.h"
#include "output_writer/python_callback/python_callback.h"
#include "output_writer/python_file/python_file.h"
#include "output_writer/paraview/paraview.h"
#include "output_writer/exfile/exfile.h"
#include "output_writer/megamol/megamol.h"

namespace OutputWriter
{

void Manager::initialize(DihuContext context, PythonConfig settings)
{
  
  LOG(TRACE) << "Entering Manager::initialize";
  outputWriter_.clear();

  //if(PyDict_Check(settings.pyObject()))
  //{
    PythonUtility::printDict(settings.pyObject());
    
    if (settings.hasKey("OutputWriter"))                    // <<--- if "settings.pyObject()" is not a Python dict, this will crash!
    {
      LOG(DEBUG) << "OutputWriter::Manager::initialize(), settings has key OutputWriter";

      // get the first value from the list
      PyObject *writerSettings = settings.getOptionListBegin<PyObject *>("OutputWriter");

      // loop over other values
      for (;
          !settings.getOptionListEnd("OutputWriter");
          settings.getOptionListNext<PyObject *>("OutputWriter", writerSettings))
      {
        if (VLOG_IS_ON(1))
        {
          VLOG(1) << "parse outputWriter, settings: " << settings.pyObject() << ", writerSettings: " << writerSettings;
        }
        PythonConfig writerConfig(settings, "OutputWriter", writerSettings);
        createOutputWriterFromSettings(context, writerConfig);
      }
    }
    else
    {
      std::vector<std::string> configKeys;
      LOG(TRACE) << "   getKeys(configKeys)";
      settings.getKeys(configKeys);
      LOG(TRACE) << "   ?";
      LOG(DEBUG) << "Config does not contain \"OutputWriter\". Keys: " << configKeys;
      //PythonUtility::printDict(settings.pyObject());
    }
    LOG(TRACE) << "Leaving Manager::initialize";
  //}
  //else
  //{
  //  LOG(ERROR) << "Can't initialize manager. PythonConfig \"" << settings << "\" is not a Python Dict.";
  //}
}

void Manager::createOutputWriterFromSettings(DihuContext context, PythonConfig settings)
{
  LOG(DEBUG) << "createOutputWriterFromSettings " << settings;
  if (settings.hasKey("format"))
  {
    // depending on type string create different OutputWriter object
    std::string typeString = settings.getOptionString("format", "none");

    LOG(DEBUG) << "add OutputWriter with format \"" << typeString << "\"";
    if (typeString == "Paraview")
    {
      outputWriter_.push_back(std::make_shared<Paraview>(context, settings));
    }
    else if (typeString == "PythonCallback")
    {
      outputWriter_.push_back(std::make_shared<PythonCallback>(context, settings));
    }
    else if (typeString == "PythonFile")
    {
      outputWriter_.push_back(std::make_shared<PythonFile>(context, settings));
    }
    else if (typeString == "Exfile" || typeString == "ExFile")
    {
      outputWriter_.push_back(std::make_shared<Exfile>(context, settings));
    }
    else if (typeString == "MegaMol")
    {
#ifdef HAVE_ADIOS
      outputWriter_.push_back(std::make_shared<MegaMol>(context, settings));
#else
      LOG(ERROR) << "Not compiled with ADIOS, but a \"MegaMol\" output writer was specified. Ignoring this output writer.";
#endif
    }
    else
    {
      LOG(WARNING) << "Unknown output writer type \"" << typeString<< "\". "
        << "Valid options are: \"Paraview\", \"PythonCallback\", \"PythonFile\", \"Exfile\", \"MegaMol\"";
    }
  }
}

bool Manager::hasOutputWriters()
{
  return !outputWriter_.empty();
}

}  // namespace
