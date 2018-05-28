#pragma once

#include <Python.h>  // has to be the first included header
#include <fstream>

#include "control/types.h"
#include "data_management/data.h"
#include "output_writer/loop_collect_mesh_names.h"

namespace OutputWriter
{

class Generic
{
public:
  //! ctor
  Generic(PyObject *specificSettings);

  //! virtual destructor to allow dynamic_pointer_cast
  virtual ~Generic();

protected:

  //! check if output should be written in this timestep and prepare filename, i.e. set filename_ from config
  template<typename DataType>
  bool prepareWrite(DataType &data, int timeStepNo = -1, double currentTime = 0.0);

  //! open file given by filename, create directory if necessary
  std::ofstream openFile(std::string filename);

  std::string filenameBase_;    ///< beginning of the file name for output file
  std::string filename_;        ///< file name with time step number
  int writeCallCount_ = 0;           ///< counter of calls to write

  int timeStepNo_;              ///< the current time step no.
  double currentTime_;          ///< the current simulation time

  PyObject *specificSettings_;    ///< the python dict containing settings relevant to this objecto
};

};  // namespace

#include "output_writer/generic.tpp"
