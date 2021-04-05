#pragma once

#include <Python.h>  // has to be the first included header

#include "cellml/source_code_generator/02_generator_openmp.h"

#include <functional>
#include <set>
#include <vc_or_std_simd.h>

#ifndef HAVE_STDSIMD      // only if we are using Vc, it is not necessary for std::simd
#include <Vc/Allocator>
#endif

class CellmlSourceCodeGeneratorVc :
  public CellmlSourceCodeGeneratorOpenMp
{
public:
  //! constructor of parent class
  using CellmlSourceCodeGeneratorOpenMp::CellmlSourceCodeGeneratorOpenMp;

  //! write the source file with explicit vectorization using Vc
  //! The file contains the source for the total solve the rhs computation
  void generateSourceFileFastMonodomain(std::string outputFilename, bool approximateExponentialFunction, 
                                        int approximateExponentialFunctionSeriesIndex);

protected:

  //! create Vc constructs for scalar functions (ternary operator) and pow/exp functions
  //! if useVc is false, no Vc:: constructs will be employed
  void preprocessCode(std::set<std::string> &helperFunctions, bool useVc = true);

  //! define "pow" and "exponential" helper functions
  std::string defineHelperFunctions(std::set<std::string> &helperFunctions, bool approximateExponentialFunction, 
                                    int approximateExponentialFunctionSeriesIndex, bool useVc, bool useReal = true);

  //! Write the source file with explicit vectorization using Vc
  //! The file contains the source for only the rhs computation
  void generateSourceFileVc(std::string outputFilename, bool approximateExponentialFunction, int approximateExponentialFunctionSeriesIndex, bool useAoVSMemoryLayout=false);

  bool preprocessingDone_ = false;      //< if preprocessing of the code tree has been done already
  std::string helperFunctionsCode_;     //< code with all helper functions like pow, exponential
};
