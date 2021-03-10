#pragma once

#include <Python.h>  // has to be the first included header

#include "cellml/source_code_generator/03_generator_vc.h"

class CellmlSourceCodeGeneratorGpu :
  public CellmlSourceCodeGeneratorVc
{
public:
  //! constructor of parent class
  using CellmlSourceCodeGeneratorVc::CellmlSourceCodeGeneratorVc;

  //! generate source code for use in the GPU fast monodomain solver
  void generateSourceFastMonodomainGpu(bool approximateExponentialFunction, int nFibersToCompute, int nInstancesToComputePerFiber,
                                       int nParametersPerInstance, bool hasAlgebraicsForTransfer,
                                       std::string &headerCode, std::string &mainCode);

protected:

  //! write the source file with openmp support
  void generateSourceFileGpu(std::string outputFilename);
};
