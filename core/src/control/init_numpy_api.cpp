// https://stackoverflow.com/questions/47026900/pyarray-check-gives-segmentation-fault-with-cython-c

#include "use_numpy.h"

//! function that initializes numpy before execution of any code
void initNumpy()
{

#ifdef HAVE_NUMPYC
  if (_import_array() < 0)
  {
    PyErr_Print();
    PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
  }
#endif

}
