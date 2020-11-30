#pragma once

#include <Python.h>  // has to be the first included header

#include "cellml/source_code_generator/00_source_code_generator_base.h"

class CellmlSourceCodeGeneratorSimd : public CellmlSourceCodeGeneratorBase
{
public:
  //! constructor of parent class
  using CellmlSourceCodeGeneratorBase::CellmlSourceCodeGeneratorBase;

protected:

  //! write the source file with openmp pragmas in struct-of-array memory ordering
  //! that will be autovectorized by the compiler
  void generateSourceFileSimd(std::string outputFilename);

};
