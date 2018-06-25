#include "solver/solver.h"

#include "utility/python_utility.h"

namespace Solver
{

Solver::Solver(PyObject *specificSettings) : specificSettings_(specificSettings)
{

}

bool Solver::configEquals(PyObject* config)
{
  if (config && specificSettings_)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;
  
    return PyObject_RichCompareBool(specificSettings_, config, Py_EQ);
  }
  return false;
}



};  // namespace