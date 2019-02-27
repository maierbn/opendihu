#include "control/dihu_context.h"

#include <Python.h>  // this has to be the first included header

void DihuContext::initializeLogging(int argc, char *argv[])
{
  START_EASYLOGGINGPP(argc, argv);
/*
  std::ifstream file("logging.conf");
  if (!file.is_open())
  {
    // if file does not exist, create it
    std::ofstream out("logging.conf");
    if (!out.is_open())
    {
      LOG(ERROR) << "Could not open logging file for output";
    }
    out << R"(
* GLOBAL:
   FORMAT               =  "INFO : %msg"
   FILENAME             =  "/tmp/logs/my.log"
   ENABLED              =  true
   TO_FILE              =  true
   TO_STANDARD_OUTPUT   =  true
   SUBSECOND_PRECISION  =  1
   PERFORMANCE_TRACKING =  false
   MAX_LOG_FILE_SIZE    =  2097152 ## 2MB - Comment starts with two hashes (##)
   LOG_FLUSH_THRESHOLD  =  100 ## Flush after every 100 logs
* DEBUG:
   FORMAT               = "DEBUG: %msg"
* WARNING:
   FORMAT               = "WARN : %loc %func: Warning: %msg"
* ERROR:
   FORMAT               = "ERROR: %loc %func: Error: %msg"
* FATAL:
   FORMAT               = "FATAL: %loc %func: Fatal error: %msg"
    )";
  }
  file.close();

  el::Configurations conf("logging.conf");
*/

// color codes: https://github.com/shiena/ansicolor/blob/master/README.md
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_LIGHT_GRAY    "\x1b[90m"
#define ANSI_COLOR_LIGHT_WHITE    "\x1b[97m"
#define ANSI_COLOR_RESET   "\x1b[0m"

  std::string separator(80, '_');
  el::Configurations conf;
  conf.setToDefault();

  // set prefix for output that includes current rank no
  std::string prefix;
  if (nRanksCommWorld_ > 1)
  {
    std::stringstream s;
    s << ownRankNoCommWorld_ << "/" << nRanksCommWorld_ << " ";
    prefix = s.str();
  }
  
#ifdef NDEBUG      // if release
  if (nRanksCommWorld_ > 1)
  {
    conf.setGlobally(el::ConfigurationType::Format, prefix+": %msg");
  }
  else
  {
    conf.setGlobally(el::ConfigurationType::Format, "%msg");
  }
#else
  conf.setGlobally(el::ConfigurationType::Format, prefix+"INFO : %msg");
#endif

  // set location of log files
  std::string logFilesPath = "/tmp/logs/";   // must end with '/'
  if (nRanksCommWorld_ > 1)
  {
    std::stringstream s;
    s << logFilesPath << ownRankNoCommWorld_ << "_opendihu.log";
    conf.setGlobally(el::ConfigurationType::Filename, s.str());

    // truncate logfile
    std::ofstream logfile(s.str().c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
    logfile.close();
  }
  else
  {
    std::string logFilename = logFilesPath+"opendihu.log";
    conf.setGlobally(el::ConfigurationType::Filename, logFilename);

    // truncate logfile
    std::ofstream logfile(logFilename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
    logfile.close();
  }

  conf.setGlobally(el::ConfigurationType::Enabled, "true");
  conf.setGlobally(el::ConfigurationType::ToFile, "true");
  conf.setGlobally(el::ConfigurationType::ToStandardOutput, "true");

  // set format of outputs
  conf.set(el::Level::Debug, el::ConfigurationType::Format, prefix+"DEBUG: %msg");
  conf.set(el::Level::Trace, el::ConfigurationType::Format, prefix+"TRACE: %msg (" ANSI_COLOR_LIGHT_GRAY "%func" ANSI_COLOR_RESET " at "
    ANSI_COLOR_LIGHT_GRAY "%loc" ANSI_COLOR_RESET ")");
  conf.set(el::Level::Verbose, el::ConfigurationType::Format, ANSI_COLOR_LIGHT_WHITE "" + prefix+"VERB%vlevel: %msg" ANSI_COLOR_RESET);
  conf.set(el::Level::Warning, el::ConfigurationType::Format,
  //         prefix+"WARN : %loc %func: \n" ANSI_COLOR_YELLOW "Warning: " ANSI_COLOR_RESET "%msg");
           prefix+ANSI_COLOR_YELLOW "Warning: " ANSI_COLOR_RESET "%msg");

  conf.set(el::Level::Error, el::ConfigurationType::Format,
           prefix+"ERROR: %loc %func: \n" ANSI_COLOR_RED "Error: %msg" ANSI_COLOR_RESET);

  conf.set(el::Level::Fatal, el::ConfigurationType::Format,
           "FATAL: %loc %func: \n"+std::string(ANSI_COLOR_MAGENTA)+prefix+separator
           +"\n\nFatal error: %msg\n"+separator+ANSI_COLOR_RESET+"\n");

  // disable output for ranks != 0
  if (ownRankNoCommWorld_ > 0)
  {
    conf.set(el::Level::Info, el::ConfigurationType::Enabled, "false");
    conf.set(el::Level::Warning, el::ConfigurationType::Enabled, "false");
  }

  //el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);

//#ifdef NDEBUG      // if release
//  conf.set(el::Level::Debug, el::ConfigurationType::Enabled, "false");
//  std::cout<< "DISABLE Debug" << std::endl;
//#endif

  // reconfigure all loggers
  el::Loggers::reconfigureAllLoggers(conf);
  el::Loggers::removeFlag(el::LoggingFlag::AllowVerboseIfModuleNotSpecified);
}
