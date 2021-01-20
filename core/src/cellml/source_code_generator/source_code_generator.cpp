#include "cellml/source_code_generator/source_code_generator.h"

#include <Python.h>  // has to be the first included header

void CellmlSourceCodeGenerator::generateSourceFile(std::string outputFilename, std::string optimizationType, bool approximateExponentialFunction, int maximumNumberOfThreads)
{
  if (optimizationType == "vc")
  {
    generateSourceFileVc(outputFilename, approximateExponentialFunction);
  }
  else if (optimizationType == "simd")
  {
    generateSourceFileSimd(outputFilename);
  }
  else if (optimizationType == "openmp")
  {
    generateSourceFileOpenMP(outputFilename, maximumNumberOfThreads);
  }
  else if (optimizationType == "gpu")
  {
    generateSourceFileGpu(outputFilename);
  }
}
