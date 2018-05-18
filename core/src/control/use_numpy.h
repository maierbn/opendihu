// https://stackoverflow.com/questions/47026900/pyarray-check-gives-segmentation-fault-with-cython-c

#ifdef HAVE_NUMPYC  // if the build system has successfully built the numpy c-API
  #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

  #include <Python.h>  // has to be the first included header

  //now, everything is setup, just include the numpy-arrays:
  #include <numpy/arrayobject.h>

#endif  // HAVE_NUMPYC

//! function that initializes numpy before execution of any code
void initNumpy();
