#include "control/dihu_context.h"

#include <Python.h>  // this has to be the first included header
#include <python_home.h>  // defines PYTHON_HOME_DIRECTORY
#include "control/diagnostic_tool/performance_measurement.h"
#include "utility/python_capture_stderr.h"

void DihuContext::initializePython(int argc, char *argv[], bool explicitConfigFileGiven)
{
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

  // add the emb module that captures stderr to the existing table of built-in modules
  PyImport_AppendInittab("emb", emb::PyInit_emb);

  // initialize python
  Py_Initialize();

  PyEval_InitThreads();

  //VLOG(4) << "PyEval_ReleaseLock()";
  //PyEval_ReleaseLock();

  Py_SetStandardStreamEncoding(NULL, NULL);

  // import emb module which captures stderr
  PyImport_ImportModule("emb");

  // get standard python path
  wchar_t *standardPythonPathWChar = Py_GetPath();
  std::wstring standardPythonPath(standardPythonPathWChar);

  VLOG(1) << "standard python path: " << standardPythonPath;

  // set python path
  //std::stringstream pythonPath;
  //pythonPath << ".:" << PYTHON_HOME_DIRECTORY << "/lib/python3.6:" << PYTHON_HOME_DIRECTORY << "/lib/python3.6/site-packages:"
  //pythonPath << OPENDIHU_HOME << "/scripts:" << OPENDIHU_HOME << "/scripts/geometry_manipulation";
  //VLOG(1) << "python path: " << pythonPath.str();
  //const wchar_t *pythonPathWChar = Py_DecodeLocale(pythonPath.str().c_str(), NULL);
  //Py_SetPath((wchar_t *)pythonPathWChar);

  // adjust PYTHONPATH
  std::stringstream codeForPythonPath;
  codeForPythonPath << "import sys" << std::endl
    << "sys.path.append('" << OPENDIHU_HOME << "/scripts" << "')" << std::endl
    << "sys.path.append('" << OPENDIHU_HOME << "/scripts/geometry_manipulation" << "')" << std::endl;
  PyRun_SimpleString(codeForPythonPath.str().c_str());


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
  // set own rank no and number of ranks in parameters
  Control::PerformanceMeasurement::setParameter("rankNo", ownRankNoCommWorld_);
  Control::PerformanceMeasurement::setParameter("nRanks", nRanksCommWorld_);

  // convert to wchar_t
  std::stringstream rankNoStr, nRanksStr;
  rankNoStr << ownRankNoCommWorld_;
  nRanksStr << nRanksCommWorld_;
  argumentToConfigWChar[nArgumentsToConfig-2] = Py_DecodeLocale(rankNoStr.str().c_str(), NULL);
  argumentToConfigWChar[nArgumentsToConfig-1] = Py_DecodeLocale(nRanksStr.str().c_str(), NULL);

  if (VLOG_IS_ON(1) && pythonConfig_.pyObject())
  {
    PythonUtility::printDict(pythonConfig_.pyObject());
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

}

bool DihuContext::loadPythonScriptFromFile(std::string filename)
{
  // initialize python interpreter

  std::ifstream file(filename);
  if (!file.is_open())
  {
    LOG(WARNING) << "Could not open settings file \"" <<filename << "\".";
    return false;
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

    std::string pythonScriptFilePath;

    // get current working directory
    char currentWorkingDirectory[PATH_MAX+1];
    if (getcwd(currentWorkingDirectory, sizeof(currentWorkingDirectory)))
    {

      // set a command that assigns the absolute path to the current settings file to the __file__ attribute
      std::stringstream commandToSetFile;
      // if filename of settings is an absolute path, starting with '/'
      if (filename[0] == '/')
      {
        commandToSetFile << "__file__ = '" << filename << "' "<< std::endl;
      }
      else
      {
        commandToSetFile << "__file__ = '" << currentWorkingDirectory << "/" << filename << "' "<< std::endl;
      }
      LOG(DEBUG) << "run " << commandToSetFile.str();

      int ret = PyRun_SimpleString(commandToSetFile.str().c_str());
      PythonUtility::checkForError();

      // if there was an error in the python code
      if (ret != 0)
      {
        PyErr_Print();
      }
    }
    else
    {
      LOG(WARNING) << "Could not set __file__ constant.";
    }

    loadPythonScript(fileContents);
  }
  return true;
}

void DihuContext::loadPythonScript(std::string text)
{
  pythonScriptText_ = text;
  LOG(TRACE) << "loadPythonScript(" << pythonScriptText_.substr(0,std::min(std::size_t(80),pythonScriptText_.length())) << ")";

  // execute python code
  int ret = 0;
  std::string errorBuffer;
  std::string stdoutBuffer;
  LOG(INFO) << std::string(40, '-') << " begin python output " << std::string(40, '-');
  try
  {
    // check if numpy module could be loaded
    PyObject *numpyModule = PyImport_ImportModule("numpy");
    if (numpyModule == NULL)
    {
      // get standard python path
      wchar_t *standardPythonPathWChar = Py_GetPath();
      std::wstring standardPythonPath(standardPythonPathWChar);
      LOG(ERROR) << "Failed to import numpy. \n Python home directory: \"" << PYTHON_HOME_DIRECTORY
        << "\", Standard python path: " << standardPythonPath;

      wchar_t *homeWChar = Py_GetPythonHome();
      char *home = Py_EncodeLocale(homeWChar, NULL);
      LOG(ERROR) << "python home: " << home;

      wchar_t *pathWChar = Py_GetPath();
      char *path = Py_EncodeLocale(pathWChar, NULL);
      LOG(ERROR) << "python path: " << path;

      const char *version = Py_GetVersion();
      LOG(ERROR) << "python version: " << version;
    }

    // add callback function to capture stderr buffer
    emb::stderr_write_type stderrWrite = [&errorBuffer] (std::string s) {errorBuffer += s; };
    emb::set_stderr(stderrWrite);

    // add callback function to capture stdout buffer
    emb::stdout_write_type stdoutWrite = [&stdoutBuffer] (std::string s) {std::cout << s; stdoutBuffer += s; };
    emb::set_stdout(stdoutWrite);

    // execute config script
    ret = PyRun_SimpleString(pythonScriptText_.c_str());

    emb::reset_stderr();
    emb::reset_stdout();

    LOG(DEBUG) << stdoutBuffer;

    PythonUtility::checkForError();
  }
  catch(...)
  {
  }
  LOG(INFO) << std::string(40, '-') << "- end python output -" << std::string(40, '-');

  // if there was an error in the python code
  if (ret != 0)
  {
    if (PyErr_Occurred())
    {
      // print error message and exit
      PyErr_Print();
      LOG(FATAL) << "An error occured in the python config.\n" << errorBuffer;
    }

    PyErr_Print();
    LOG(FATAL) << "An error occured in the python config.\n" << errorBuffer;
  } else if (!errorBuffer.empty()) {
    LOG(WARNING) << "The python config wrote to stderr.\n";
    LOG(INFO) << std::string(37, '-') << " begin python error output " << std::string(37, '-');
    LOG(INFO) << errorBuffer;
    LOG(INFO) << std::string(37, '-') << "- end python error output -" << std::string(37, '-');
  }

  // load main module and extract config
  PyObject *mainModule = PyImport_AddModule("__main__");
  PyObject *config = PyObject_GetAttrString(mainModule, "config");
  VLOG(4) << "create pythonConfig_ (initialize ref to 1)";


  // check if type is valid
  if (config == NULL || !PyDict_Check(config))
  {
    LOG(FATAL) << "Python config file does not contain a dict named \"config\".";
  }

  pythonConfig_.setPyObject(config);

  parseGlobalParameters();
}

void DihuContext::parseGlobalParameters()
{
  // parse scenario name
  std::string scenarioName = pythonConfig_.getOptionString("scenarioName", "");
  Control::PerformanceMeasurement::setParameter("scenarioName", scenarioName);

  // parse desired log file format
  std::string logFormat = pythonConfig_.getOptionString("logFormat", "csv");
  if (logFormat == "csv")
  {
    setLogFormat(logFormatCsv);
  }
  else if (logFormat == "json")
  {
    setLogFormat(logFormatJson);
  }
  else
  {
    LOG(ERROR) << "Unknown option for \"logFormat\": \"" << logFormat << "\". Use one of \"csv\" or \"json\". Falling back to \"csv\".";
    setLogFormat(logFormatCsv);
  }

  // parse all keys under meta and add forward them directly to the log file
  // These parameters are not used by opendihu but can hold information that the
  // user wants to have in the log file.
  if (pythonConfig_.hasKey("meta"))
  {

    // loop over entries of python dict "meta"
    std::string keyString("meta");
    std::pair<std::string,std::string> dictItem
      = pythonConfig_.getOptionDictBegin<std::string,std::string>(keyString);

    for (; !pythonConfig_.getOptionDictEnd(keyString);
        pythonConfig_.getOptionDictNext<std::string,std::string>(keyString, dictItem))
    {
      // the key is the name of the field
      std::string key = dictItem.first;

      // the value is the payload string
      std::string value = dictItem.second;

      std::stringstream logKey;
      logKey << "meta_" << key;
      LOG(DEBUG) << "key[" << key << "] value [" << value << "]";
      Control::PerformanceMeasurement::setParameter(logKey.str(), value);
    }
  }

  // filename for solver structure diagram
  solverStructureDiagramFile_ = pythonConfig_.getOptionString("solverStructureDiagramFile", "solver_structure.txt");
  if (solverStructureDiagramFile_ == "None")
    solverStructureDiagramFile_ = "";
}
