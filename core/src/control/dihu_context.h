#pragma once

#include <petscsys.h>
#include <Python.h>
#include <list>
#include <memory>

#include <output_writer/generic.h>
#include <data_management/data.h>

class DihuContext
{
public:
  ///! constructor, initialize context, parse command line parameters and input file
  DihuContext(int argc, char *argv[]);
  
  ///! constructor for test cases
  DihuContext(int argc, char *argv[], std::string pythonSettings);
  
  ///! return the top-level python config object
  PyObject *getPythonConfig();
  
  ///! get reference to a PetscErrorCode temporary variable to be used to assign petsc error codes
  PetscErrorCode &ierr();
  
  ///! call all output writers to write output, timeStepNo of -1 means no time step number in output filename
  void writeOutput(Data::Data &problemData, PyObject *specificSettings, 
                   int timeStepNo = -1, double currentTime = 0.0);
  
  ///! destructor
  ~DihuContext();
  
private:
  ///! read in file and execute python script and store global variables
  void loadPythonScriptFromFile(std::string filename);
  
  ///! execute python script and store global variables
  void loadPythonScript(std::string text);
  
  ///! initiaize the easylogging++ framework
  void initializeLogging();
  
  ///! create output writer from settings
  void initializeOutputWriter();
  
  ///! helper function that creates on outputWriter
  void createOutputWriterFromSettings(PyObject *dict);
  
  PyObject *pythonConfig_;    ///< the top level python config dictionary
  
  PetscErrorCode ierr_;     ///< temporary variable for petsc error codes
  std::list<std::unique_ptr<OutputWriter::Generic>> outputWriter_;    ///< list of active output writers
};