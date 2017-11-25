#include "control/dihu_context.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <list>

#include "control/python_utility.h"
#include "output_writer/paraview.h"
#include "output_writer/python.h"
#include "output_writer/callback.h"
#include "mesh/mesh_manager.h"

#include "Python.h"
#include "easylogging++.h"

//INITIALIZE_EASYLOGGINGPP

std::list<std::unique_ptr<OutputWriter::Generic>> DihuContext::outputWriter_;
std::shared_ptr<MeshManager> DihuContext::meshManager_ = nullptr;
bool DihuContext::initialized_ = false;
 
DihuContext::DihuContext(int argc, char *argv[]) :
  pythonConfig_(NULL)
{
  LOG(TRACE) << "DihuContext constructor";

  if (!initialized_)
  {
    // load configuration from file if it exits
    initializeLogging(argc, argv);
    
    // initialize MPI, this is necessary to be able to call PetscFinalize without MPI shutting down
    MPI_Init(&argc, &argv);
    
    // initialize PETSc
    PetscInitialize(&argc, &argv, NULL, "This is an opendihu application.");
    
    // determine settings filename
    std::string filename = "settings.py";
    
    if (argc > 1)
    {
      if (argv[1][0] != '-')    // do not consider command line arguments starting with '-' as input settings file
        filename = argv[1];
    }
    
    
    LOG(TRACE) << "initialize python";
    
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
    int ret = system("printf %s $(find /usr -name \"libpython*.a\" | head -n 1 | sed 's/\\/lib.*//') > tmp");
    if (ret == 0)
    {
      std::ifstream f("tmp");
      std::stringstream s;
      s << f.rdbuf();
      std::remove("tmp");
      
      const char *pythonSearchPath = s.str().c_str();
      LOG(DEBUG) << "Set python search path to \""<<pythonSearchPath<<"\".";
      
      Py_SetPythonHome((char *)pythonSearchPath);
    }
    
    char const *programName = "dihu";
    Py_SetProgramName((char *)programName);  /* optional but recommended */
    
    Py_Initialize();
    
    char *home = Py_GetPythonHome();
    LOG(DEBUG) << "Python home: " << home;
    
    loadPythonScriptFromFile(filename);
    initializeOutputWriter();
    
    initialized_ = true;
  }

  if (!meshManager_)
  {
    LOG(TRACE) << "create meshManager_";
    meshManager_ = std::make_shared<MeshManager>(*this);
  }
}  

DihuContext::DihuContext(int argc, char *argv[], std::string pythonSettings) : DihuContext(argc, argv)
{
  loadPythonScript(pythonSettings);
}

PyObject* DihuContext::getPythonConfig() const
{
  //if (!pythonConfig_)
  //  LOG(FATAL) << "Python config is not available!";
  return pythonConfig_;
}

std::shared_ptr<MeshManager> DihuContext::meshManager() const
{
  return meshManager_;
}

DihuContext DihuContext::operator[](std::string keyString) const
{
  int argc = 0;
  char **argv = NULL;
  DihuContext dihuContext(argc, argv);
  if (PythonUtility::containsKey(pythonConfig_, keyString))
  {
    dihuContext.pythonConfig_ = PythonUtility::getOptionPyObject(pythonConfig_, keyString);
    Py_XINCREF(dihuContext.pythonConfig_);
  }
  else
  {
    dihuContext.pythonConfig_ = pythonConfig_;
    Py_XINCREF(dihuContext.pythonConfig_);
    LOG(WARNING) << "Dict does not contain key \""<<keyString<<"\".";
  }
  LOG(TRACE) << "DihuContext::operator[](\""<<keyString<<"\")";
  
  dihuContext.initializeOutputWriter();
  
  return dihuContext;
}

void DihuContext::loadPythonScriptFromFile(std::string filename)
{
  LOG(TRACE)<<"loadPythonScriptFromFile";
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
  LOG(TRACE)<<"loadPythonScript";
  
  // execute python code
  int ret = 0;
  LOG(INFO)<<std::string(80, '-');
  try
  {
    ret = PyRun_SimpleString(text.c_str());
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

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
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
}

void DihuContext::initializeOutputWriter()
{
  if (PythonUtility::containsKey(pythonConfig_, "OutputWriter"))
  {
    // get the first value from the list
    PyObject *writerSettings = PythonUtility::getOptionListBegin<PyObject *>(pythonConfig_, "OutputWriter");
  
    // loop over other values
    for (;
        !PythonUtility::getOptionListEnd(pythonConfig_, "OutputWriter");
        PythonUtility::getOptionListNext<PyObject *>(pythonConfig_, "OutputWriter", writerSettings))
    {
      createOutputWriterFromSettings(writerSettings);
    }
  }
#if 0
  PyObject *specificSettings = PythonUtility::getOptionPyObject(topLevelSettings, "OutputWriter");
  
  if (specificSettings)
  {
    if (PyList_Check(specificSettings))
    {
      // type is list
      for(int i=0; i<PyList_Size(specificSettings); i++)
      {
        // extract entry
        PyObject *listEntry = PyList_GetItem(specificSettings, (Py_ssize_t)i);
        
        if (!PyDict_Check(listEntry))
        {
          LOG(WARNING) << "Discard non-dict entry of list in OutputWriter";
          continue;
        }
        createOutputWriterFromSettings(listEntry);
      }
    }
    else if(PyDict_Check(specificSettings))
    {
      createOutputWriterFromSettings(specificSettings);
    }
    else
    {
      LOG(WARNING) << "Entry \"OutputWriter\" in config has to be a list of dicts or a dict.";
    }
  }
  else
  {
    LOG(WARNING) << "config does not contain key \"OutputWriter\"!";
  }
#endif
}

void DihuContext::createOutputWriterFromSettings(PyObject *dict)
{
  PyObject *key = PyString_FromString("format");
  if (PyDict_Contains(dict, key))
  {
    PyObject *type = PyDict_GetItem(dict, key);
    if (PyString_Check(type))
    {
      std::string typeString = PyString_AsString(type);
      if (typeString == "Paraview")
      {
        outputWriter_.push_back(std::make_unique<OutputWriter::Paraview>(dict));
      }
      else if(typeString == "Python")
      {
        outputWriter_.push_back(std::make_unique<OutputWriter::Python>(dict));
      }
      else if(typeString == "Callback")
      {
        outputWriter_.push_back(std::make_unique<OutputWriter::Callback>(dict));
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

DihuContext::~DihuContext()
{
  // do not finalize Python because otherwise tests keep crashing
  Py_CLEAR(pythonConfig_);
  //Py_Finalize();

  // do not finalize Petsc because otherwise there can't be multiple DihuContext objects for testing
  //PetscErrorCode ierr;
  //ierr = PetscFinalize(); CHKERRV(ierr);
  //MPI_Finalize();
}

void DihuContext::writeOutput(Data::Data &problemData, int timeStepNo, double currentTime) const
{
  for(auto &outputWriter : outputWriter_)
  {
    outputWriter->write(problemData, timeStepNo, currentTime);
  }
}