#include "control/dihu_context.h"

#include <Python.h>  // this has to be the first included header
#include <python_home.h>  // defines PYTHON_HOME_DIRECTORY

#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <list>

#include "utility/python_utility.h"
#include "output_writer/paraview/paraview.h"
#include "output_writer/python_callback/python_callback.h"
#include "output_writer/python_file/python_file.h"
#include "output_writer/exfile/exfile.h"
#include "mesh/mesh_manager.h"
#include "solver/solver_manager.h"

#include "easylogging++.h"
#include "control/use_numpy.h"

//INITIALIZE_EASYLOGGINGPP

std::shared_ptr<Mesh::Manager> DihuContext::meshManager_ = nullptr;
std::shared_ptr<Solver::Manager> DihuContext::solverManager_ = nullptr;

bool DihuContext::initialized_ = false;

// copy-constructor
DihuContext::DihuContext(const DihuContext &rhs)
{
  pythonConfig_ = rhs.pythonConfig_;
  VLOG(4) << "Py_XINCREF(pythonConfig_)";
  Py_XINCREF(pythonConfig_);
}

DihuContext::DihuContext(int argc, char *argv[], bool settingsFromFile) :
  pythonConfig_(NULL)
{
  LOG(TRACE) << "DihuContext constructor";

  if (!initialized_)
  {
    // load configuration from file if it exits
    initializeLogging(argc, argv);

    // initialize MPI, this is necessary to be able to call PetscFinalize without MPI shutting down
    MPI_Init(&argc, &argv);

    // configure PETSc to abort on errorm
    PetscOptionsSetValue(NULL, "-on_error_abort", "");

    // initialize PETSc
    PetscInitialize(&argc, &argv, NULL, "This is an opendihu application.");

    // determine settings filename
    std::string filename = "settings.py";

    // check if the first command line argument is *.py, only then it is treated as config file
    bool explicitConfigFileGiven = false;
    if (argc > 1)
    {
      std::string firstArgument = argv[1];
      if (firstArgument.rfind(".py") == firstArgument.size() - 3)
      {
        explicitConfigFileGiven = true;
        filename = argv[1];
      }
      else
      {
        LOG(INFO) << "First command line argument does not have suffix *.py, not considering it as config file!";
      }
    }

    LOG(TRACE) << "initialize python";
#if 1

    // set program name of python script
    char const *programName = "dihu";
    VLOG(4) << "Py_SetProgramName(" << std::string(programName) << ")";

#if PY_MAJOR_VERSION >= 3
    wchar_t *programNameWChar = Py_DecodeLocale(programName, NULL);
    Py_SetProgramName(programNameWChar);  /* optional but recommended */
#else
    Py_SetProgramName((char *)programName);  /* optional but recommended */
#endif

#if PY_MAJOR_VERSION >= 3
    // set python home and path, apparently this is not needed
#if 1
    VLOG(1) << "python home directory: \"" << PYTHON_HOME_DIRECTORY << "\"";
    std::string pythonSearchPath = PYTHON_HOME_DIRECTORY;
    //std::string pythonSearchPath = std::string("/store/software/opendihu/dependencies/python/install");
    const wchar_t *pythonSearchPathWChar = Py_DecodeLocale(pythonSearchPath.c_str(), NULL);
    Py_SetPythonHome((wchar_t *)pythonSearchPathWChar);

    //std::string pythonPath = ".:" PYTHON_HOME_DIRECTORY "/lib/python3.6:" PYTHON_HOME_DIRECTORY "/lib/python3.6/site-packages";
    //const wchar_t *pythonPathWChar = Py_DecodeLocale(pythonPath.c_str(), NULL);
    //Py_SetPath((wchar_t *)pythonPathWChar);
#endif

#else   // python 2.7

    // set python home and path

    // find location where libpython2.7.a is located and extract the beginning of the path until '/lib' is encountered
    // this serves as the PYTHONHOME value that will be set via Py_SetPythonHome.
    // possible examples:
    //  /usr/lib/x86_64-linux-gnu   -> PYTHONHOME=/usr
    //  /usr/lib/python2.7/config-x86_64-linux-gnu -> PYTHONHOME=/usr
    // explanation of the command:
    //   printf %s $(..): print without newline
    //   find /usr -name "libpython*.a": search for libpython*.a under the /usr directory
    //   head -n 1: take first line
    //   sed 's/\/lib.*//': remove everything after "/lib" is found
    //   > tmp: write to file "tmp"
    //int ret = system("printf %s $(find /usr -name \"libpython*.a\" | head -n 1 | sed 's/\\/lib.*//') > tmp");

    int ret = 0;
    if (ret == 0)
    {
      std::ifstream f("tmp");
      std::stringstream s;
      s << f.rdbuf();
      std::remove("tmp");

      std::string pythonHome = s.str();
      if (pythonHome.empty())
        pythonHome = "/usr";

      const char *pythonSearchPath = pythonHome.c_str();
      LOG(DEBUG) << "Set python search path to \""<<pythonSearchPath<<"\".";

      VLOG(4) << "Py_SetPythonHome(" << pythonHome << ")";
    }
    Py_SetPythonHome((char *)pythonSearchPath);
#endif

    // initialize python
    VLOG(4) << "Py_Initialize()";
    Py_Initialize();

    VLOG(4) << "Py_SetStandardStreamEncoding()";
    Py_SetStandardStreamEncoding(NULL, NULL);

    // pass on command line arguments to python config script

    // determine if the first argument (argv[1]) is *.py, then it is also discarded
    // always remove the first argument, which is the name of the executable
    int numberArgumentsToRemove = (explicitConfigFileGiven? 2: 1);

    int argcReduced = argc - numberArgumentsToRemove;
    VLOG(4) << "argcReduced: " << argcReduced << ", numberArgumentsToRemove: " << numberArgumentsToRemove;

    char **argvReduced = new char *[argcReduced];
#if PY_MAJOR_VERSION >= 3
    wchar_t **argvReducedWChar = new wchar_t *[argcReduced];
#endif

    for (int i=0; i<argcReduced; i++)
    {
      argvReduced[i] = argv[i+numberArgumentsToRemove];
      LOG(DEBUG) << "argv[" << i << "]=" << argvReduced[i];
#if PY_MAJOR_VERSION >= 3
      argvReducedWChar[i] = Py_DecodeLocale(argvReduced[i], NULL);
#endif
    }

    // pass reduced list of command line arguments to python script
#if PY_MAJOR_VERSION >= 3
    PySys_SetArgvEx(argcReduced, argvReducedWChar, 0);
#else
    PySys_SetArgvEx(argcReduced, argvReduced, 0);
#endif



#if PY_MAJOR_VERSION >= 3
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
#else
    // check python home directory for debugging output
    char *home = Py_GetPythonHome();
    VLOG(2) << "Python home: " << home;
#endif

    // load python script
    if(settingsFromFile)
    {
      loadPythonScriptFromFile(filename);
    }

#endif
    initialized_ = true;
  }
  // if this is the first constructed DihuContext object, create global objects mesh manager and solver manager
  if (!meshManager_)
  {
    VLOG(2) << "create meshManager_";
    meshManager_ = std::make_shared<Mesh::Manager>(pythonConfig_);
  }

  if (!solverManager_)
  {
    VLOG(2) << "create solverManager_";
    solverManager_ = std::make_shared<Solver::Manager>(pythonConfig_);
  }
}

DihuContext::DihuContext(int argc, char *argv[], std::string pythonSettings) : DihuContext(argc, argv, false)
{
  loadPythonScript(pythonSettings);
  PythonUtility::printDict(pythonConfig_);

  VLOG(2) << "recreate meshManager";
  meshManager_ = nullptr;
  meshManager_ = std::make_shared<Mesh::Manager>(pythonConfig_);
  solverManager_ = nullptr;
  solverManager_ = std::make_shared<Solver::Manager>(pythonConfig_);

}

PyObject* DihuContext::getPythonConfig() const
{
  //if (!pythonConfig_)
  //  LOG(FATAL) << "Python config is not available!";
  return pythonConfig_;
}

std::shared_ptr<Mesh::Manager> DihuContext::meshManager() const
{
  return meshManager_;
}

std::shared_ptr<Solver::Manager> DihuContext::solverManager() const
{
  return solverManager_;
}

DihuContext DihuContext::operator[](std::string keyString)
{
  int argc = 0;
  char **argv = NULL;
  DihuContext dihuContext(argc, argv);

  // if requested child context exists in config
  if (PythonUtility::hasKey(pythonConfig_, keyString))
  {
    dihuContext.pythonConfig_ = PythonUtility::getOptionPyObject(pythonConfig_, keyString);
    VLOG(4) << "Py_XINCREF(dihuContext.pythonConfig_)";
    Py_XINCREF(dihuContext.pythonConfig_);
  }
  else
  {
    // if config does not contain the requested child dict, create the needed context from the same level in config
    dihuContext.pythonConfig_ = pythonConfig_;
    VLOG(4) << "Py_XINCREF(dihuContext.pythonConfig_)";
    Py_XINCREF(dihuContext.pythonConfig_);
    LOG(WARNING) << "Dict does not contain key \""<<keyString<<"\".";
  }
  LOG(TRACE) << "DihuContext::operator[](\""<<keyString<<"\")";

  return dihuContext;
}

DihuContext DihuContext::createSubContext(PyObject *settings)
{
  int argc = 0;
  char **argv = NULL;
  DihuContext dihuContext(argc, argv);

  dihuContext.pythonConfig_ = settings;
  VLOG(4) << "Py_XINCREF(dihuContext.pythonConfig_)";
  Py_XINCREF(dihuContext.pythonConfig_);

  return dihuContext;
}


void DihuContext::loadPythonScriptFromFile(std::string filename)
{
  // initialize python interpreter

  std::ifstream file(filename);
  if (!file.is_open())
  {
    LOG(ERROR)<<"Could not open settings file \""<<filename<<"\".";
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

    LOG(INFO)<<"File \""<<filename<<"\" loaded.";

    loadPythonScript(fileContents);
  }
}

void DihuContext::loadPythonScript(std::string text)
{
  LOG(TRACE)<<"loadPythonScript("<<text.substr(0,std::min(std::size_t(80),text.length()))<<")";

  // execute python code
  int ret = 0;
  LOG(INFO)<<std::string(80, '-');
  try
  {
    LOG(DEBUG) << "run python script";
    PyObject *numpyModule = PyImport_ImportModule("numpy");
    if (numpyModule == NULL)
      LOG(DEBUG) << "failed to import numpy";
    ret = PyRun_SimpleString(text.c_str());
    if (PyErr_Occurred())
    {
      PyErr_Print();
    }
  }
  catch(...)
  {
  }
  LOG(INFO)<<std::string(80, '-');

  // if there was an error in the python code
  if (ret != 0)
  {
    if (PyErr_Occurred())
    {
      // print error message and exit
      PyErr_Print();
      exit(0);
    }
    exit(0);
  }

  // load main module
  PyObject *mainModule = PyImport_AddModule("__main__");
  pythonConfig_ = PyObject_GetAttrString(mainModule, "config");
  VLOG(4) << "create pythonConfig_ (initialize ref to 1)";


  // check if type is valid
  if (pythonConfig_ == NULL || !PyDict_Check(pythonConfig_))
  {
    LOG(ERROR)<<"Python config file does not contain a dict named \"config\".";
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

  conf.setGlobally(el::ConfigurationType::Format, "INFO : %msg");
  conf.setGlobally(el::ConfigurationType::Filename, "/tmp/logs/my.log");
  conf.setGlobally(el::ConfigurationType::Enabled, "true");
  conf.setGlobally(el::ConfigurationType::ToFile, "true");
  conf.setGlobally(el::ConfigurationType::ToStandardOutput, "true");

  // set format of outputs
  conf.set(el::Level::Debug, el::ConfigurationType::Format, "DEBUG: %msg");
  conf.set(el::Level::Trace, el::ConfigurationType::Format, "TRACE: %msg");
  conf.set(el::Level::Verbose, el::ConfigurationType::Format, ANSI_COLOR_LIGHT_WHITE "VERB%vlevel: %msg" ANSI_COLOR_RESET);
  conf.set(el::Level::Warning, el::ConfigurationType::Format,
           "WARN : %loc %func: \n" ANSI_COLOR_YELLOW "Warning: " ANSI_COLOR_RESET "%msg");

  conf.set(el::Level::Error, el::ConfigurationType::Format,
           "ERROR: %loc %func: \n" ANSI_COLOR_RED "Error: %msg" ANSI_COLOR_RESET);

  conf.set(el::Level::Fatal, el::ConfigurationType::Format,
           std::string(ANSI_COLOR_MAGENTA)+"FATAL: %loc %func: \n"+separator
           +"\nFatal error: %msg\n"+separator+ANSI_COLOR_RESET+"\n");

  //el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);

//#ifdef NDEBUG      // if release
//  conf.set(el::Level::Debug, el::ConfigurationType::Enabled, "false");
//  std::cout<<"DISABLE Debug"<<std::endl;
//#endif

  // reconfigure all loggers
  el::Loggers::reconfigureAllLoggers(conf);
  el::Loggers::removeFlag(el::LoggingFlag::AllowVerboseIfModuleNotSpecified);
}

DihuContext::~DihuContext()
{
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
