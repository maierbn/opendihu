#include <iostream>
#include <cstdlib>

#include "opendihu.h"
#include "utility/svd_utility.h"

int main(int argc, char *argv[])
{
  //SvdUtility::getSVD();
  SvdUtility::readCSV("trajectories.csv");
  return EXIT_SUCCESS;
}