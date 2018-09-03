#pragma once

#include <Python.h>  // this has to be the first included header
#include <iostream>
#include <vector>

//! assert that the file given by filename has exactly the content given in referenceContent or referenceContent2, fail the test otherwise
void assertFileMatchesContent(std::string filename, std::string referenceContent, std::string referenceContent2="-");

//! check that the parallel and serial output files contain the same data, using the script "validate_parallel.py"
void assertParallelEqualsSerialOutputFiles(std::vector<std::string> &outputFilesToCheck);

extern int nFails;
