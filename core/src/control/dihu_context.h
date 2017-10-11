#pragma once

#include <petscsys.h>

class DihuContext
{
public:
  /// constructor, initialize context, parse command line parameters and input file
  DihuContext(int argc, char *argv[]);
  
  ///get reference to a PetscErrorCode temporary variable to be used to assign petsc error codes
  PetscErrorCode &ierr();
  
private:
  
  /// execute python script and store global variables
  void loadPythonScript(std::string filename);
  
  PetscErrorCode ierr_;     ///< temporary variable for petsc error codes
};