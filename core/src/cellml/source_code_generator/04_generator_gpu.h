#pragma once

#include <Python.h>  // has to be the first included header

#include "cellml/source_code_generator/03_generator_vc.h"

class CellmlSourceCodeGeneratorGpu :
  public CellmlSourceCodeGeneratorVc
{
public:
  //! constructor of parent class
  using CellmlSourceCodeGeneratorVc::CellmlSourceCodeGeneratorVc;

  // to be defined
};
