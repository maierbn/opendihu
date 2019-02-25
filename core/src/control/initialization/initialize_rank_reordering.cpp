#include "control/dihu_context.h"

#include <Python.h>  // this has to be the first included header
#include <python_home.h>  // defines PYTHON_HOME_DIRECTORY
#include "control/performance_measurement.h"

void DihuContext::initializeRankReordering(int argc, char *argv[])
{
  rankReorderingEnabled_ = false;

  // check if the --rank_reordering flag is given on the command line

  for (int i = 0; i < argc; i++)
  {
    std::string argument(argv[i]);
    if (argument == "--rank_reordering")
    {
      rankReorderingEnabled_ = true;
      break;
    }
  }

  if (!rankReorderingEnabled_)
    return;

  // parse subdomains from arguments
  for (int i = 0; i < argc; i++)
  {
    std::string argument(argv[i]);
    if (argument == "--n_subdomains")
    {
      nSubdomainsForRankReordering_[0] = atoi(argv[i+1]);
      nSubdomainsForRankReordering_[1] = atoi(argv[i+2]);
      nSubdomainsForRankReordering_[2] = atoi(argv[i+3]);
      LOG(DEBUG) << "parsed nSubdomainsForRankReordering_: " << nSubdomainsForRankReordering_;
      break;
    }
  }
}

void DihuContext::reorderRankNoCommWorld(int &rankNo)
{
  if (rankReorderingEnabled_)
  {
    // determine x,y,z coordinates
    int subdomainCoordinateX = rankNo % nSubdomainsForRankReordering_[0];
    int subdomainCoordinateY = (int)(rankNo / nSubdomainsForRankReordering_[0]) % nSubdomainsForRankReordering_[1];
    int subdomainCoordinateZ = (int)(rankNo / (nSubdomainsForRankReordering_[0]*nSubdomainsForRankReordering_[1]));

    int subdomainCoordinateXNew = subdomainCoordinateZ;
    int subdomainCoordinateYNew = subdomainCoordinateX;
    int subdomainCoordinateZNew = subdomainCoordinateY;

    // assemble new rank number
    rankNo = subdomainCoordinateZNew*nSubdomainsForRankReordering_[1]*nSubdomainsForRankReordering_[0]
      + subdomainCoordinateYNew*nSubdomainsForRankReordering_[0] + subdomainCoordinateXNew;
  }
}
