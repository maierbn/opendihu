#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // Given a binary file with fibers scale all points and write output file
  if (argc != 4)
  {
    std::cout << "usage: scale <input_filename> <output_filename> <scaling_factor>";
    return EXIT_SUCCESS;
  }

  std::string inputFilename = argv[1];
  std::string outputFilename = argv[2];
  double scalingFactor = atof(argv[3]);

  Postprocessing::scaleFibersInFile(inputFilename, outputFilename, scalingFactor);

  return EXIT_SUCCESS;
}
