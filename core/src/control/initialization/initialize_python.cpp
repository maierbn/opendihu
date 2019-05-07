#include "control/dihu_context.h"

#include <Python.h>  // this has to be the first included header
#include <control/python_home.h>  // defines PYTHON_HOME_DIRECTORY
#include "control/performance_measurement.h"
//#include <dlfcn.h>

void DihuContext::initializePython(int argc, char *argv[], bool explicitConfigFileGiven)
{
  LOG(TRACE) << "initialize python with explicitConfigFileGiven=" << explicitConfigFileGiven;

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

//#ifdef __PGI
//  //workaround (see https://stackoverflow.com/questions/11842920/undefined-symbol-pyexc-importerror-when-embedding-python-in-c)
//  void*const libpythonHandle = dlopen("/afs/.mathematik.uni-stuttgart.de/home/cmcs/share/environment-modules/Packages/python/python-3.6.5/lib/libpython3.so", RTLD_LAZY | RTLD_GLOBAL);
//  // dlclose counterpart in dihu_context.cpp NOT implemented.   
//#endif
  
  // initialize python
  Py_Initialize();

  PyEval_InitThreads();

  //VLOG(4) << "PyEval_ReleaseLock()";
  //PyEval_ReleaseLock();

  Py_SetStandardStreamEncoding(NULL, NULL);

  // get standard python path
  wchar_t *standardPythonPathWChar = Py_GetPath();
  std::wstring standardPythonPath(standardPythonPathWChar);

  VLOG(1) << "standard python path: " << standardPythonPath;

  // set python path
  std::stringstream pythonPath;
  //pythonPath << ".:" << PYTHON_HOME_DIRECTORY << "/lib/python3.6:" << PYTHON_HOME_DIRECTORY << "/lib/python3.6/site-packages:"
  //pythonPath << OPENDIHU_HOME << "/scripts:" << OPENDIHU_HOME << "/scripts/geometry_manipulation";
  //VLOG(1) << "python path: " << pythonPath.str();
  //const wchar_t *pythonPathWChar = Py_DecodeLocale(pythonPath.str().c_str(), NULL);
  //Py_SetPath((wchar_t *)pythonPathWChar);


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
  
  LOG(DEBUG) << "try printDict:";
  PythonUtility::printDict(pythonConfig_.pyObject());

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
  LOG(TRACE) << "python version: " << version;

  const char *platform = Py_GetPlatform();
  VLOG(2) << "python platform: " << platform;

  const char *compiler = Py_GetCompiler();
  VLOG(2) << "python compiler: " << compiler;

  const char *buildInfo = Py_GetBuildInfo();
  VLOG(2) << "python buildInfo: " << buildInfo;

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
    LOG(TRACE) << "File \"" <<filename << "\" loaded.";

    loadPythonScript(fileContents);
  }
}

void DihuContext::loadPythonScript(std::string text)
{
  pythonScriptText_ = text;
  LOG(TRACE) << "loadPythonScript(" << pythonScriptText_.substr(0,std::min(std::size_t(80),pythonScriptText_.length())) << ")";

  // execute python code
  int ret = 0;
  LOG(INFO) << std::string(80, '-');
  try
  { LOG(TRACE) << "check import";
    // check if numpy module could be loaded
    //LOG(INFO) << "PyRun_SimpleString(\"import numpy\"): " << PyRun_SimpleString("import numpy");
    PyObject *numpyModule = PyImport_ImportModule("numpy");
    if (numpyModule == NULL)
    {
      LOG(ERROR) << "Failed to import numpy.";
      if (PyErr_Occurred())
      {
        LOG(INFO) << "PyErr occurred.";
        PyErr_Print();
      }
      else{
        LOG(INFO) << "PyErr_Occurred() = NULL";}
    }
    LOG(TRACE) << "execute config script";
    // execute config script
    ret = PyRun_SimpleString(pythonScriptText_.c_str());
    //LOG(INFO) << "PyRun_SimpleString(arg) with arg: " << pythonScriptText_.c_str();

    PythonUtility::checkForError();
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

  LOG(TRACE) << "load main module";
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

  // parse scenario name
  std::string scenarioName = "";
  if (pythonConfig_.hasKey("scenarioName"))
  {
    scenarioName = pythonConfig_.getOptionString("scenarioName", "");
  }
  Control::PerformanceMeasurement::setParameter("scenarioName", scenarioName);
}
