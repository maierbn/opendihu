#include "output_writer/manager.h"

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "output_writer/python_callback/python_callback.h"
#include "output_writer/python_file/python_file.h"
#include "output_writer/paraview/paraview.h"
#include "output_writer/exfile/exfile.h"
#include "output_writer/megamol/megamol.h"

namespace OutputWriter
{

void Manager::initialize(DihuContext context, PythonConfig settings, std::shared_ptr<Partition::RankSubset> rankSubset)
{
  std::vector<int> outputFileNo;
  if (!outputWriter_.empty())
  {
    for (std::list<std::shared_ptr<Generic>>::iterator iter = outputWriter_.begin(); iter != outputWriter_.end(); iter++)
    {
      outputFileNo.push_back((*iter)->outputFileNo());
    }
  }

  outputWriter_.clear();

  //VLOG(3) << "initializeOutputWriter, settings=" << settings;
  //PythonUtility::printDict(settings.pyObject());

  if (settings.hasKey("OutputWriter"))
  {
    LOG(DEBUG) << "OutputWriter::Manager::initialize(), settings has key OutputWriter";

    // get the first value from the list
    PyObject *writerSettings = settings.template getOptionListBegin<PyObject *>("OutputWriter");

    // loop over other values
    for (;
        !settings.getOptionListEnd("OutputWriter");
        settings.template getOptionListNext<PyObject *>("OutputWriter", writerSettings))
    {
      if (VLOG_IS_ON(1))
      {
        VLOG(1) << "parse outputWriter, settings: " << settings.pyObject() << ", writerSettings: " << writerSettings;
      }
      PythonConfig writerConfig(settings, "OutputWriter", writerSettings);
      createOutputWriterFromSettings(context, writerConfig, rankSubset);
    }
  }
  else if (settings.hasKey("OutputWriters"))
  {
    LOG(ERROR) << "You wrote \"OutputWriters\" but it should be \"OutputWriter\".";
  }
  else
  {
    std::vector<std::string> configKeys;
    settings.getKeys(configKeys);

#ifndef NDEBUG
    VLOG(1) << "Config does not contain \"OutputWriter\". Keys: " << configKeys;
    //PythonUtility::printDict(settings.pyObject());
#endif
  }


  // if there were outputFileNo stored from earlier output writers, restore these numbers
  if (!outputWriter_.empty() && !outputFileNo.empty())
  {
    int i = 0;
    for (std::list<std::shared_ptr<Generic>>::iterator iter = outputWriter_.begin(); iter != outputWriter_.end(); iter++)
    {
      (*iter)->setOutputFileNo(outputFileNo[i]);

      if (i < outputFileNo.size())
        i++;
    }
  }

}

void Manager::createOutputWriterFromSettings(DihuContext context, PythonConfig settings, std::shared_ptr<Partition::RankSubset> rankSubset)
{
  LOG(DEBUG) << "createOutputWriterFromSettings " << settings;
  if (settings.hasKey("format"))
  {
    // depending on type string create different OutputWriter object
    std::string typeString = settings.getOptionString("format", "none");

    LOG(DEBUG) << "add OutputWriter with format \"" << typeString << "\"";
    if (typeString == "Paraview")
    {
      outputWriter_.push_back(std::make_shared<Paraview>(context, settings, rankSubset));
    }
    else if (typeString == "PythonCallback")
    {
      outputWriter_.push_back(std::make_shared<PythonCallback>(context, settings, rankSubset));
    }
    else if (typeString == "PythonFile")
    {
      outputWriter_.push_back(std::make_shared<PythonFile>(context, settings, rankSubset));
    }
    else if (typeString == "Exfile" || typeString == "ExFile")
    {
      outputWriter_.push_back(std::make_shared<Exfile>(context, settings, rankSubset));
    }
    else if (typeString == "MegaMol")
    {
#ifdef HAVE_ADIOS
      outputWriter_.push_back(std::make_shared<MegaMol>(context, settings, rankSubset));
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

//! get the filename of the first output writer
std::string Manager::filename()
{
  if (outputWriter_.empty())
    return std::string("");

  return outputWriter_.front()->filenameBase();
}

//! set the filename for the first output writer
void Manager::setFilename(std::string filename)
{
  if (!outputWriter_.empty())
  {
    outputWriter_.front()->setFilenameBase(filename);
  }
}

}  // namespace
