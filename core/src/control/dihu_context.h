#pragma once

#include <petscsys.h>
#include <Python.h>
#include <list>
#include <memory>

#include <output_writer/generic.h>

class DihuContext
{
public:
  ///! constructor, initialize context, parse command line parameters and input file
  DihuContext(int argc, char *argv[]);
  
  ///! return the top-level python config object
  PyObject *getPythonConfig();
  
  ///! get reference to a PetscErrorCode temporary variable to be used to assign petsc error codes
  PetscErrorCode &ierr();
  
  ///! call all output writers to write output 
  void writeOutput(ProblemData &problemData, PyObject *specificSettings);
  
  
  ///! destructor
  ~DihuContext();
  
private:
  
  ///! execute python script and store global variables
  void loadPythonScript(std::string filename);
  
  ///! create output writer from settings
  void initializeOutputWriter();
  
  ///! helper function that creates on outputWriter
  void createOutputWriterFromSettings(PyObject *dict);
  
  PyObject *pythonConfig_;    ///< the top level python config dictionary
  
  PetscErrorCode ierr_;     ///< temporary variable for petsc error codes
  std::list<std::unique_ptr<OutputWriter::Generic>> outputWriter_;    ///< list of active output writers
};