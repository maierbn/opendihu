#pragma once

#include <petscsys.h>
#include <Python.h>

class DihuContext
{
public:
  /// constructor, initialize context, parse command line parameters and input file
  DihuContext(int argc, char *argv[]);
  
  /// return the top-level python config object
  PyObject *getPythonConfig();
  
  ///get reference to a PetscErrorCode temporary variable to be used to assign petsc error codes
  PetscErrorCode &ierr();
  
  /// destructor
  ~DihuContext();
  
private:
  
  /// execute python script and store global variables
  void loadPythonScript(std::string filename);
  
  PyObject *pythonConfig_;    ///< the top level python config dictionary
  
  PetscErrorCode ierr_;     ///< temporary variable for petsc error codes
};