#include "control/dihu_context.h"

#include <Python.h>  // this has to be the first included header

#ifdef HAVE_MEGAMOL
#include "Console.h"     // contains megamol_main
#endif


#include <adios2.h>

void DihuContext::initializeAdios(int argc, char *argv[])
{
  LOG(DEBUG) << "initializeAdios";
  adios_ = std::make_shared<adios2::ADIOS>(MPI_COMM_WORLD);
  io_ = std::make_shared<adios2::IO>(adios_->DeclareIO("Output"));
  assert(io_);
}

void DihuContext::initializeMegaMol(int argc, char *argv[])
{
  // extract MegaMol arguments from config
  if (pythonConfig_.hasKey("MegaMolArguments"))
  {
    initializeAdios(argc, argv);
#ifdef HAVE_MEGAMOL


    std::string megamolArgumentsOption = pythonConfig_.getOptionString("MegaMolArguments", "");

    megamolArguments_.clear();
    // convert options from a string to a vector of strings
    size_t pos = 0;
    while (pos < megamolArgumentsOption.length())
    {
      size_t newPos = megamolArgumentsOption.find(" ", pos);
      if (newPos == std::string::npos)
      {
        megamolArguments_.push_back(megamolArgumentsOption.substr(pos));
        break;
      }

      megamolArguments_.push_back(megamolArgumentsOption.substr(pos, newPos-pos));
      pos = megamolArgumentsOption.find(" ", pos+1)+1;
    }

    // prepare arguments
    int megamolArgc = megamolArguments_.size()+1;
    megamolArgv_.resize(megamolArgc);
    static std::string commandName = "mmconsole";
    megamolArgv_[0] = (char *)commandName.c_str();

    for (int i = 0; i < megamolArguments_.size(); i++)
    {
      megamolArgv_[i+1] = (char *)megamolArguments_[i].c_str();
    }

    LOG(DEBUG) << "start MegaMol with arguments " << megamolArguments_;

    //start megamol main method
    megamolThread_ = std::make_shared<std::thread>(megamol_main, megamolArgc, megamolArgv_.data());
#else
    LOG(ERROR) << "Not compiled with MegaMol, but \"MegaMolArguments\" are given.";
#endif
  }
}
