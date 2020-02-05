#pragma once

#include <Python.h>  // has to be the first included header

#include "cellml/source_code_generator/04_generator_gpu.h"

class CellmlSourceCodeGenerator :
  public CellmlSourceCodeGeneratorGpu
{
public:
  //! constructor
  using CellmlSourceCodeGeneratorGpu::CellmlSourceCodeGeneratorGpu;

  //! generate the source file according to optimizationType
  //! Possible values are: simd vc openmp
  //! @param approximateExponentialFunction If the exp()-Function should be approximated by the n=1024th series term
  void generateSourceFile(std::string outputFilename, std::string optimizationType, bool approximateExponentialFunction);
};
