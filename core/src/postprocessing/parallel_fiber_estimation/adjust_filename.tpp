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

  LOG(INFO) << "adjust filename [" << filename << "], nFibersX: " << nFibersX << ", path: [" << path << "], fileBase: [" << fileBase << "]";

  // create directory and wait until system has created it
  if (path != "")
  {
    int ret = system((std::string("mkdir -p ")+path).c_str());
    if (ret != 0)
      LOG(WARNING) << "Creation of directory \"" <<path<< "\" failed.";
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
  }

  // check if filename contains x
  bool xFound = false;
  bool filenameHasXFormat = false;
  int suffixPos = 0;  // character position of the suffix
  int xNumbersBeginPos = 0;
  int firstDigitPos = -1;

  for(int i = 0; i < fileBase.size(); i++)
  {
    if (isdigit(fileBase[i]))
    {
      if (firstDigitPos == -1)
        firstDigitPos = i;
    }
    else
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
        xNumbersBeginPos = firstDigitPos;
      }
      firstDigitPos = -1;
    }
  }

  LOG(DEBUG) << " filenameHasXFormat: " << filenameHasXFormat;

  if (filenameHasXFormat)
  {
    std::stringstream newFilename;
    newFilename << path << fileBase.substr(0,xNumbersBeginPos) << nFibersX << "x" << nFibersX  << fileBase.substr(suffixPos);
    filename = newFilename.str();
    LOG(INFO) << "newFilename: [" << filename << "]";
    return true;
  }
  else
  {
    return false;
  }

}

} // namespace
