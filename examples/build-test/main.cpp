#include <iostream>
#include <cstdlib>

#include "control/simulation.h"

int main(int argc, char *argv[])
{
  std::cout << "hu" << std::endl;
  Simulation simulation;
  simulation.debug();
  
  return EXIT_SUCCESS;
}