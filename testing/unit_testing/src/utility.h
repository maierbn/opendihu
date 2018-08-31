#pragma once

#include <Python.h>  // this has to be the first included header
#include <iostream>

//! assert that the file given by filename has exactly the content given in referenceContent or referenceContent2, fail the test otherwise
void assertFileMatchesContent(std::string filename, std::string referenceContent, std::string referenceContent2="-");

extern int nFails;
