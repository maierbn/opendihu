#pragma once

#include <Python.h>

#include "control/types.h"
#include "data_management/data.h"

namespace OutputWriter
{
 
class Generic
{
public:
  Generic(PyObject *specificSettings);
  
  //! write output file, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = 0.0);
protected:
 
  std::string filenameBase_;    ///< beginning of the file name for output file
  std::string filename_;        ///< file name with time step number
  int writeCallCount_ = 0;           ///< counter of calls to write
  
  PyObject *specificSettings_;    ///< the python dict containing settings relevant to this objecto
};

};  // namespace 

#include "output_writer/generic.tpp"