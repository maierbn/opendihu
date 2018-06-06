#pragma once

#include <Python.h>  // has to be the first included header

#include <Python.h>
#include <vector>

#include "control/dihu_context.h"
#include "output_writer/manager.h"
#include "cellml/00_cellml_adapter_base.h"

/** This is a class that handles the right hand side routine of CellML
 * It needs access on the system size and some C-code of the right hand side.
 * It provides a method to call the right hand side.
 */
template <int nStates>
class RhsRoutineHandler:
  public CellmlAdapterBase<nStates>
{
public:
  
  //! constructor
  using CellmlAdapterBase<nStates>::CellmlAdapterBase;
  
protected:
  //! given a normal cellml source file for rhs routine create a second file for multiple instances. @return: if successful
  bool createSimdSourceFile(std::string &simdSourceFilename);

  //! initialize the rhs routine, either directly from a library or compile it
  void initializeRhsRoutine();

  void (*rhsRoutine_)(double, double *, double *, double *, double *);   ///< function pointer to a rhs function that is passed as dynamic library, computes rates and intermediate values from states. The parameters are: VOI, STATES, RATES, WANTED, KNOWN, (VOI: unclear, what it means)
  void (*rhsRoutineSimd_)(void *context, double *, double *, double *, double *);  ///< same functionality as rhsRoutine, however, can compute several instances of the problem in parallel. Data is assumed to contain values for a state contiguously, e.g. (state[1], state[1], state[1], state[2], state[2], state[2], ...). The first parameter is a this pointer

};

#include "cellml/01_rhs_routine_handler.tpp"