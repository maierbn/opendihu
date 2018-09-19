#include "control/dihu_context.h"

#include <Python.h>  // this has to be the first included header
#include <python_home.h>  // defines PYTHON_HOME_DIRECTORY
#include <omp.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <list>
#include <petscvec.h>

#include "utility/python_utility.h"
#include "output_writer/paraview/paraview.h"
#include "output_writer/python_callback/python_callback.h"
#include "output_writer/python_file/python_file.h"
#include "output_writer/exfile/exfile.h"
#include "mesh/mesh_manager.h"
#include "solver/solver_manager.h"
#include "partition/partition_manager.h"

#include "easylogging++.h"
#include "control/use_numpy.h"
#include "utility/mpi_utility.h"

//INITIALIZE_EASYLOGGINGPP

std::shared_ptr<Mesh::Manager> DihuContext::meshManager_ = nullptr;
//std::shared_ptr<Solver::Manager> DihuContext::solverManager_ = nullptr;
std::map<int, std::shared_ptr<Solver::Manager>> DihuContext::solverManagerForThread_;
std::shared_ptr<Partition::Manager> DihuContext::partitionManager_ = nullptr;

bool DihuContext::initialized_ = false;
int DihuContext::nObjects_ = 0;   ///< number of objects of DihuContext, if the last object gets destroyed, call MPI_Finalize

// copy-constructor
DihuContext::DihuContext(const DihuContext &rhs)
{
  nObjects_++;
  doNotFinalizeMpi_ = rhs.doNotFinalizeMpi_;
  pythonConfig_ = rhs.pythonConfig_;
  Py_XINCREF(pythonConfig_);
}

DihuContext::DihuContext(int argc, char *argv[], bool doNotFinalizeMpi, bool settingsFromFile) :
  pythonConfig_(NULL), doNotFinalizeMpi_(doNotFinalizeMpi)
{
  nObjects_++;
  LOG(TRACE) << "DihuContext constructor";

  if (!initialized_)
  {
    // initialize MPI, this is necessary to be able to call PetscFinalize without MPI shutting down
    MPI_Init(&argc, &argv);

    // load configuration from file if it exits
    initializeLogging(argc, argv);

    // configure PETSc to abort on errorm
    PetscOptionsSetValue(NULL, "-on_error_abort", "");

    // initialize PETSc
    PetscInitialize(&argc, &argv, NULL, "This is an opendihu application.");

    // parallel debugging barrier
    bool enableDebuggingBarrier = false;
    PetscErrorCode ierr = PetscOptionsHasName(NULL, NULL, "-pause", (PetscBool *)&enableDebuggingBarrier); CHKERRV(ierr);

    if (enableDebuggingBarrier)
    {
      MPIUtility::gdbParallelDebuggingBarrier();
    }

    // determine settings filename
    std::string filename = "settings.py";

    // check if the first command line argument is *.py, only then it is treated as config file
    bool explicitConfigFileGiven = false;
    if (argc > 1 && settingsFromFile)
    {
      std::string firstArgument = argv[1];
      if (firstArgument.rfind(".py") == firstArgument.size() - 3)
      {
        explicitConfigFileGiven = true;
        filename = argv[1];
      }
      else
      {
        LOG(ERROR) << "First command line argument does not have suffix *.py, not considering it as config file!";
      }
    }

    LOG(TRACE) << "initialize python";

    // set program name of python script
    char const *programName = "opendihu";

    wchar_t *programNameWChar = Py_DecodeLocale(programName, NULL);
    Py_SetProgramName(programNameWChar);  /* optional but recommended */

    // set python home and path, apparently this is not needed
    VLOG(1) << "python home directory: \"" << PYTHON_HOME_DIRECTORY << "\"";
    std::string pythonSearchPath = PYTHON_HOME_DIRECTORY;
    //std::string pythonSearchPath = std::string("/store/software/opendihu/dependencies/python/install");
    const wchar_t *pythonSearchPathWChar = Py_DecodeLocale(pythonSearchPath.c_str(), NULL);
    Py_SetPythonHome((wchar_t *)pythonSearchPathWChar);

    //std::string pythonPath = ".:" PYTHON_HOME_DIRECTORY "/lib/python3.6:" PYTHON_HOME_DIRECTORY "/lib/python3.6/site-packages";
    //const wchar_t *pythonPathWChar = Py_DecodeLocale(pythonPath.c_str(), NULL);
    //Py_SetPath((wchar_t *)pythonPathWChar);

    // initialize python
    Py_Initialize();

    PyEval_InitThreads();
    
    //VLOG(4) << "PyEval_ReleaseLock()";
    //PyEval_ReleaseLock();

    Py_SetStandardStreamEncoding(NULL, NULL);

    // pass on command line arguments to python config script

    // determine if the first argument (argv[1]) is *.py, then it is also discarded
    // always remove the first argument, which is the name of the executable
    int numberArgumentsToRemove = (explicitConfigFileGiven? 2: 1);

    // add the own rank no and the number of ranks at the end as command line arguments

    int nArgumentsToConfig = argc - numberArgumentsToRemove + 2;
    VLOG(4) << "nArgumentsToConfig: " << nArgumentsToConfig << ", numberArgumentsToRemove: " << numberArgumentsToRemove;

    char **argvReduced = new char *[nArgumentsToConfig];
    wchar_t **argumentToConfigWChar = new wchar_t *[nArgumentsToConfig];

    // set given command line arguments
    for (int i=0; i<nArgumentsToConfig-2; i++)
    {
      argvReduced[i] = argv[i+numberArgumentsToRemove];
      argumentToConfigWChar[i] = Py_DecodeLocale(argvReduced[i], NULL);
    }

    // add rank no and nRanks
    // get own rank no and number of ranks
    int rankNo, nRanks;
    MPIUtility::handleReturnValue (MPI_Comm_rank(MPI_COMM_WORLD, &rankNo));
    MPIUtility::handleReturnValue (MPI_Comm_size(MPI_COMM_WORLD, &nRanks));

    // convert to wchar_t
    std::stringstream rankNoStr, nRanksStr;
    rankNoStr << rankNo;
    nRanksStr << nRanks;
    argumentToConfigWChar[nArgumentsToConfig-2] = Py_DecodeLocale(rankNoStr.str().c_str(), NULL);
    argumentToConfigWChar[nArgumentsToConfig-1] = Py_DecodeLocale(nRanksStr.str().c_str(), NULL);

    if (VLOG_IS_ON(1) && pythonConfig_)
    {
      PythonUtility::printDict(pythonConfig_);
    }

    // pass reduced list of command line arguments to python script
    PySys_SetArgvEx(nArgumentsToConfig, argumentToConfigWChar, 0);

    // check different python setting for debugging
    wchar_t *homeWChar = Py_GetPythonHome();
    char *home = Py_EncodeLocale(homeWChar, NULL);
    VLOG(2) << "python home: " << home;

    wchar_t *pathWChar = Py_GetPath();
    char *path = Py_EncodeLocale(pathWChar, NULL);
    VLOG(2) << "python path: " << path;

    wchar_t *prefixWChar = Py_GetPrefix();
    char *prefix = Py_EncodeLocale(prefixWChar, NULL);
    VLOG(2) << "python prefix: " << prefix;

    wchar_t *execPrefixWChar = Py_GetExecPrefix();
    char *execPrefix = Py_EncodeLocale(execPrefixWChar, NULL);
    VLOG(2) << "python execPrefix: " << execPrefix;

    wchar_t *programFullPathWChar = Py_GetProgramFullPath();
    char *programFullPath = Py_EncodeLocale(programFullPathWChar, NULL);
    VLOG(2) << "python programFullPath: " << programFullPath;

    const char *version = Py_GetVersion();
    VLOG(2) << "python version: " << version;

    const char *platform = Py_GetPlatform();
    VLOG(2) << "python platform: " << platform;

    const char *compiler = Py_GetCompiler();
    VLOG(2) << "python compiler: " << compiler;

    const char *buildInfo = Py_GetBuildInfo();
    VLOG(2) << "python buildInfo: " << buildInfo;

    // load python script
    if(settingsFromFile)
    {
      loadPythonScriptFromFile(filename);
    }

    initialized_ = true;
  }
  
  // if this is the first constructed DihuContext object, create global objects partition manager, mesh manager and solver manager
  if (!partitionManager_)
  {
    VLOG(2) << "create partitionManager_";
    partitionManager_ = std::make_shared<Partition::Manager>(pythonConfig_);
  }

  if (!meshManager_)
  {
    VLOG(2) << "create meshManager_";
    meshManager_ = std::make_shared<Mesh::Manager>(pythonConfig_);
    meshManager_->setPartitionManager(partitionManager_);
  }
  
  if (solverManagerForThread_.empty())
  {
    VLOG(2) << "create solverManagerForThread_";
    // create solver manager for thread 0
    solverManagerForThread_[0] = std::make_shared<Solver::Manager>(pythonConfig_);
    solverManagerForThread_[1] = std::make_shared<Solver::Manager>(pythonConfig_);
  }
}

DihuContext::DihuContext(int argc, char *argv[], std::string pythonSettings, bool doNotFinalizeMpi) :
  DihuContext(argc, argv, doNotFinalizeMpi, false)
{
  nObjects_++;
  // This constructor is called when creating the context object from unit tests.

  // load python config script as given by parameter
  loadPythonScript(pythonSettings);
  if (VLOG_IS_ON(1))
  {
    PythonUtility::printDict(pythonConfig_);
  }

  partitionManager_ = nullptr;
  partitionManager_ = std::make_shared<Partition::Manager>(pythonConfig_);
  
  VLOG(2) << "recreate meshManager";
  meshManager_ = nullptr;
  meshManager_ = std::make_shared<Mesh::Manager>(pythonConfig_);
  meshManager_->setPartitionManager(partitionManager_);
  
  // create solver manager for thread 0
  solverManagerForThread_.clear();
  solverManagerForThread_[0] = std::make_shared<Solver::Manager>(pythonConfig_);
  
}

PyObject* DihuContext::getPythonConfig() const
{
  //if (!pythonConfig_)
  //  LOG(FATAL) << "Python config is not available!";
  return pythonConfig_;
}

int DihuContext::ownRankNo()
{
  int rankNo;
  MPIUtility::handleReturnValue (MPI_Comm_rank(MPI_COMM_WORLD, &rankNo));
  return rankNo;
}

std::shared_ptr<Mesh::Manager> DihuContext::meshManager() const
{
  return meshManager_;
}

std::shared_ptr<Partition::Manager> DihuContext::partitionManager() const
{
  return partitionManager_;
}

std::shared_ptr<Solver::Manager> DihuContext::solverManager() const
{
  // get number of omp threads
  //int nThreads = omp_get_num_threads();
  int threadId = omp_get_thread_num();
  
  if (solverManagerForThread_.find(threadId) == solverManagerForThread_.end())
  {
    LOG(DEBUG) << "create solver manager for thread " << threadId;
    // create solver manager
    solverManagerForThread_[threadId] = std::make_shared<Solver::Manager>(pythonConfig_);
    
    LOG(DEBUG) << "(done)";
  }
  else 
  {
    LOG(DEBUG) << "solver manager for thread " << threadId << " exists";
  }
  
  return solverManagerForThread_[threadId];
}

DihuContext DihuContext::operator[](std::string keyString)
{
  int argc = 0;
  char **argv = NULL;
  DihuContext dihuContext(argc, argv, doNotFinalizeMpi_);

  // if requested child context exists in config
  if (PythonUtility::hasKey(pythonConfig_, keyString))
  {
    dihuContext.pythonConfig_ = PythonUtility::getOptionPyObject(pythonConfig_, keyString);
    Py_XINCREF(dihuContext.pythonConfig_);
  }
  else
  {
    // if config does not contain the requested child dict, create the needed context from the same level in config
    dihuContext.pythonConfig_ = pythonConfig_;
    Py_XINCREF(dihuContext.pythonConfig_);
    LOG(FATAL) << "Dict does not contain key \"" <<keyString<< "\".";
  }
  LOG(TRACE) << "DihuContext::operator[](\"" <<keyString<< "\")";

  return dihuContext;
}

DihuContext DihuContext::createSubContext(PyObject *settings)
{
  int argc = 0;
  char **argv = NULL;
  DihuContext dihuContext(argc, argv);

  dihuContext.pythonConfig_ = settings;
  Py_XINCREF(dihuContext.pythonConfig_);

  return dihuContext;
}

void DihuContext::loadPythonScriptFromFile(std::string filename)
{
  // initialize python interpreter

  std::ifstream file(filename);
  if (!file.is_open())
  {
    LOG(FATAL) << "Could not open settings file \"" <<filename << "\".";
  }
  else
  {
    // reserve memory of size of file
    file.seekg(0, std::ios::end);
    size_t fileSize = file.tellg();
    std::string fileContents(fileSize, ' ');

    // reset file pointer
    file.seekg(0, std::ios::beg);

    // read in file contents
    file.read(&fileContents[0], fileSize);

    LOG(INFO) << "File \"" <<filename << "\" loaded.";

    loadPythonScript(fileContents);
  }
}

void DihuContext::loadPythonScript(std::string text)
{
  LOG(TRACE) << "loadPythonScript(" << text.substr(0,std::min(std::size_t(80),text.length())) << ")";

  // execute python code
  int ret = 0;
  LOG(INFO) << std::string(80, '-');
  try
  {
    // check if numpy module could be loaded
    PyObject *numpyModule = PyImport_ImportModule("numpy");
    if (numpyModule == NULL)
    {
      LOG(ERROR) << "Failed to import numpy.";
    }

    // execute config script
    ret = PyRun_SimpleString(text.c_str());

    if (PyErr_Occurred())
    {
      PyErr_Print();
    }
  }
  catch(...)
  {
  }
  LOG(INFO) << std::string(80, '-');

  // if there was an error in the python code
  if (ret != 0)
  {
    if (PyErr_Occurred())
    {
      // print error message and exit
      PyErr_Print();
      LOG(FATAL) << "An error occured in the python config.";
    }

    PyErr_Print();
    LOG(FATAL) << "An error occured in the python config.";
  }

  // load main module and extract config
  PyObject *mainModule = PyImport_AddModule("__main__");
  pythonConfig_ = PyObject_GetAttrString(mainModule, "config");
  VLOG(4) << "create pythonConfig_ (initialize ref to 1)";


  // check if type is valid
  if (pythonConfig_ == NULL || !PyDict_Check(pythonConfig_))
  {
    LOG(ERROR) << "Python config file does not contain a dict named \"config\".";
  }
}

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

  int nRanks;
  int rankNo;
  MPIUtility::handleReturnValue (MPI_Comm_size(MPI_COMM_WORLD, &nRanks));
  MPIUtility::handleReturnValue (MPI_Comm_rank(MPI_COMM_WORLD, &rankNo));
  
  // set prefix for output that includes current rank no
  std::string prefix;
  if (nRanks > 1)
  {
    std::stringstream s;
    s << rankNo << "/" << nRanks << " ";
    prefix = s.str();
  }
  
  conf.setGlobally(el::ConfigurationType::Format, prefix+"INFO : %msg");

  // set location of log files
  std::string logFilesPath = "/tmp/logs/";   // must end with '/'
  if (nRanks > 1)
  {
    std::stringstream s;
    s << logFilesPath << rankNo << "_opendihu.log";
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
  conf.set(el::Level::Trace, el::ConfigurationType::Format, prefix+"TRACE: %msg");
  conf.set(el::Level::Verbose, el::ConfigurationType::Format, ANSI_COLOR_LIGHT_WHITE "" + prefix+"VERB%vlevel: %msg" ANSI_COLOR_RESET);
  conf.set(el::Level::Warning, el::ConfigurationType::Format,
           prefix+"WARN : %loc %func: \n" ANSI_COLOR_YELLOW "Warning: " ANSI_COLOR_RESET "%msg");

  conf.set(el::Level::Error, el::ConfigurationType::Format,
           prefix+"ERROR: %loc %func: \n" ANSI_COLOR_RED "Error: %msg" ANSI_COLOR_RESET);

  conf.set(el::Level::Fatal, el::ConfigurationType::Format,
           "FATAL: %loc %func: \n"+std::string(ANSI_COLOR_MAGENTA)+prefix+separator
           +"\n\nFatal error: %msg\n"+separator+ANSI_COLOR_RESET+"\n");

  // disable output for ranks != 0
  if (rankNo > 0)
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

DihuContext::~DihuContext()
{
  nObjects_--;

  VLOG(1) << "~DihuContext, nObjects = " << nObjects_;
  if (nObjects_ == 1)
  {

    // After a call to MPI_Finalize we cannot call MPI_Initialize() anymore.
    // This is only a problem when the code is tested with the GoogleTest framework, because then we want to run multiple tests in one executable.
    // In this case, do not finalize MPI, but call MPI_Barrier instead which also syncs the ranks.

    if (doNotFinalizeMpi_)
    {
      LOG(INFO) << "MPI_Barrier";
      MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
      LOG(INFO) << "MPI_Finalize";
      MPI_Finalize();
    }
  }
  // do not clear pythonConfig_ here, because it crashes
  //VLOG(4) << "PY_CLEAR(PYTHONCONFIG_)";  // note: calling VLOG in a destructor is critical and can segfault
  //Py_CLEAR(pythonConfig_);



  // do not finalize Python because otherwise tests keep crashing
  //Py_Finalize();
#if PY_MAJOR_VERSION >= 3
  //Py_Finalize();
#endif

  // do not finalize Petsc because otherwise there can't be multiple DihuContext objects for testing
  //PetscErrorCode ierr;
  //ierr = PetscFinalize(); CHKERRV(ierr);
  //MPI_Finalize();
}
