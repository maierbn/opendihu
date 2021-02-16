#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>

#include "easylogging++.h"

namespace Postprocessing
{

  //! scale all points in an input bin file by the given scaling factor, create new output file
  void scaleFibersInFile(std::string inputFilename, std::string outputFilename, double scalingFactor);

} // namespace
