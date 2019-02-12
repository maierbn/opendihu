#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // Fix the invalid fibers in a given file.

  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  Postprocessing::ParallelFiberEstimation<
    BasisFunction::LagrangeOfOrder<1>
  >
  problem(settings);

  // read in the file under resultFilename and fix all invalid fibers by interpolating from neighbouring fibers
  problem.fixInvalidFibersInFile();

  return EXIT_SUCCESS;
}
