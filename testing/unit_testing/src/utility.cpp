#include "utility.h"

#include <Python.h>  // this has to be the first included header
#include <iostream>
#include <fstream>
#include "gtest/gtest.h"

void assertFileMatchesContent(std::string filename, std::string referenceContent, std::string referenceContent2)
{
  // read in generated exnode file 
  std::ifstream file(filename);

  ASSERT_TRUE(file.is_open()) << "could not open output file"; 
  
  // reserve memory of size of file
  file.seekg(0, std::ios::end);   
  size_t fileSize = file.tellg();
  std::string fileContents(fileSize, ' ');
  
  // reset file pointer
  file.seekg(0, std::ios::beg);
  
  // read in file contents
  file.read(&fileContents[0], fileSize);
  
  if (fileContents != referenceContent && fileContents != referenceContent2)
  {
    ASSERT_EQ(fileContents, referenceContent) << "file content of file \"" << filename << "\" is different (referenceContent)";
    ASSERT_EQ(fileContents, referenceContent2) << "file content of file \"" << filename << "\" is different (referenceContent2)";
  }
}
