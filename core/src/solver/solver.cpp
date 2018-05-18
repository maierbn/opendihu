#include "solver/solver.h"

namespace Solver
{

Solver::Solver(PyObject *specificSettings) : specificSettings_(specificSettings)
{

}

bool Solver::configEquals(PyObject* config)
{
  if (config && specificSettings_)
  {
    return PyObject_RichCompareBool(specificSettings_, config, Py_EQ);
  }
  return false;
}



};  // namespace