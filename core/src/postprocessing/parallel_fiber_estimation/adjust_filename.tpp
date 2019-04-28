#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

#include <sstream>

namespace Postprocessing
{

template<typename BasisFunctionType>
bool ParallelFiberEstimation<BasisFunctionType>::
adjustFilename(std::string &filename, int nFibersX)
{
  // split path and filename base
  std::string fileBase = filename;
  std::string path = "";

  if (fileBase.find("/") != std::string::npos)
  {
    path = fileBase.substr(0, fileBase.rfind("/")+1);
    fileBase = fileBase.substr(fileBase.rfind("/")+1);
  }

  // check if filename contains x
  bool xFound = false;
  bool filenameHasXFormat = false;
  int suffixPos = 0;
  for(int i = 0; i < fileBase.size(); i++)
  {
    if (!isdigit(fileBase[i]))
    {
      if (xFound)
      {
        suffixPos = i;
        filenameHasXFormat = true;
        break;
      }
      if (fileBase[i] == 'x')
      {
        xFound = true;
      }
      else
      {
        break;
      }
    }
  }

  if (filenameHasXFormat)
  {
    std::stringstream newFilename;
    newFilename << path << nFibersX << "x" << nFibersX  << fileBase.substr(suffixPos);
    filename = newFilename.str();
    return true;
  }
  else
  {
    return false;
  }

}

} // namespace
