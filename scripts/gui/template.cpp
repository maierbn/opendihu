#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
 problem(settings);
  
  // run problem
  problem.run();
  
  return EXIT_SUCCESS;
}
