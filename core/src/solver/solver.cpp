#include "solver/solver.h"

#include "utility/python_utility.h"

namespace Solver
{

Solver::Solver(PythonConfig specificSettings, std::string name) :
  specificSettings_(specificSettings), name_(name)
{
  durationLogKey_ = std::string("durationSolve_") + name_;
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

std::string Solver::name()
{
  return name_;
}

} // namespace
