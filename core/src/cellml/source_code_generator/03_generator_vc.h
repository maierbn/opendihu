#pragma once

#include <Python.h>  // has to be the first included header

#include "cellml/source_code_generator/02_generator_openmp.h"

#include <functional>
#include <Vc/Allocator>
#include <set>

class CellmlSourceCodeGeneratorVc :
  public CellmlSourceCodeGeneratorOpenMp
{
public:
  //! constructor of parent class
  using CellmlSourceCodeGeneratorOpenMp::CellmlSourceCodeGeneratorOpenMp;

  //! write the source file with explicit vectorization using Vc
  //! The file contains the source for the total solve the rhs computation
  void generateSourceFileVcFastMonodomain(std::string outputFilename, bool approximateExponentialFunction);

protected:

  //! create Vc constructs for scalar functions (ternary operator) and pow/exp functions
  void preprocessCode(std::set<std::string> &helperFunctions);

  //! define "pow" and "exponential" helper functions
  std::string defineHelperFunctions(std::set<std::string> &helperFunctions, bool approximateExponentialFunction);

  //! Write the source file with explicit vectorization using Vc
  //! The file contains the source for only the rhs computation
  void generateSourceFileVc(std::string outputFilename, bool approximateExponentialFunction);


};
