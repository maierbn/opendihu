#include "solver/solver.h"

#include "utility/python_utility.h"

namespace Solver
{

Solver::Solver(PythonConfig specificSettings, std::string name) :
  specificSettings_(specificSettings), name_(name)
{

}

bool Solver::configEquals(PythonConfig config)
{
  if (config.pyObject() != nullptr && specificSettings_.pyObject() != nullptr)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;
  
    return PyObject_RichCompareBool(specificSettings_.pyObject(), config.pyObject(), Py_EQ);
  }
  return false;
}



};  // namespace
