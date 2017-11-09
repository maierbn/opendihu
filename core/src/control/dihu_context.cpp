#include "control/dihu_context.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>

#include "control/python_utility.h"
#include "output_writer/paraview.h"
#include "output_writer/python.h"

#include "Python.h"
#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

DihuContext::DihuContext(int argc, char *argv[]) : pythonConfig_(NULL)
{
  START_EASYLOGGINGPP(argc, argv);

  // load configuration from file
  el::Configurations conf("logging.conf");
  
  // reconfigure all loggers
  el::Loggers::reconfigureAllLoggers(conf);
  
  LOG(DEBUG) << "DihuContext constructor" << std::endl;
  
  // initialize PETSc
  PetscInitialize(&argc, &argv, NULL, "This is an opendihu application.");
  
  // determine settings filename
  std::string filename = "settings.py";
  
  if (argc > 1)
  {
    filename = argv[1];
  }
  
  loadPythonScriptFromFile(filename);
}  

DihuContext::DihuContext(int argc, char *argv[], std::string pythonSettings)
{ 
  DihuContext(argc, argv);
  loadPythonScript(pythonSettings);
}

PyObject* DihuContext::getPythonConfig()
{
  //if (!pythonConfig_)
  //  LOG(FATAL) << "Python config is not available!";
  return pythonConfig_;
}

void DihuContext::loadPythonScriptFromFile(std::string filename)
{
  LOG(DEBUG)<<"loadPythonScriptFromFile";
  // initialize python interpreter
  
  char const *programName = "dihu";
  Py_SetProgramName((char *)programName);  /* optional but recommended */
  Py_Initialize();
  
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
  LOG(DEBUG)<<"loadPythonScript: "<<text;
  
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
  pythonConfig_ = PyObject_GetAttrString(mainModule,"config");
  
  // check if type is valid
  if (pythonConfig_ == NULL || !PyDict_Check(pythonConfig_))
  {
    LOG(ERROR)<<"python config file does not contain a dict named \"config\".";
  }
  else 
  {
    initializeOutputWriter();
  }
}

void DihuContext::initializeOutputWriter()
{
  PyObject *topLevelSettings = pythonConfig_;
  PyObject *specificSettings = PythonUtility::extractDict(topLevelSettings, "OutputWriter");
  
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
}

DihuContext::~DihuContext()
{
  Py_Finalize();
  if (pythonConfig_)
    Py_DECREF(pythonConfig_);

  PetscErrorCode ierr;
  ierr = PetscFinalize(); CHKERRV(ierr);
}

PetscErrorCode &DihuContext::ierr()
{
  return ierr_;
}

void DihuContext::writeOutput(Data::Data &problemData, PyObject *specificSettings, int timeStepNo, double currentTime)
{
  for(auto &outputWriter : outputWriter_)
  {
    outputWriter->write(problemData, specificSettings, timeStepNo, currentTime);
  }
}