#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"

namespace OutputWriter
{

class MegaMol : public Generic
{
public:

  //! constructor
  MegaMol(DihuContext context, PythonConfig specificSettings);

  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);

};

} // namespace

#include "output_writer/megamol/megamol.tpp"
