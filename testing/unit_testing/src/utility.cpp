#include "utility.h"

#include <Python.h>  // this has to be the first included header
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <thread>

#include "gtest/gtest.h"
#include "easylogging++.h"

int nFails = 0;   // global number of test failures

double parseNumber(std::string::iterator &iterFileContents, std::string::iterator iterFileContentsEnd)
{
  if (iterFileContents == iterFileContentsEnd)
    return 0;
  
  std::string::iterator iterFileContentsBegin = iterFileContents;

  // parse number in file Contents
  std::string numberFileContents;
  while ((isdigit(*iterFileContents) || *iterFileContents == '.' || *iterFileContents == '-' || *iterFileContents == 'e') && iterFileContents != iterFileContentsEnd)
  {
    numberFileContents += *iterFileContents;
    iterFileContents++;
  }
  
  if (numberFileContents.empty() || numberFileContents == ".")
    return 0;
  
  double number = 0;
  try 
  {
    std::stod(numberFileContents, nullptr);
  }
  catch (std::exception &e)
  {
    std::string s;
    for(std::string::iterator iterFileContents2 = iterFileContentsBegin; iterFileContents2 != iterFileContentsEnd; iterFileContents2++)
    {
      s += *iterFileContents2;
    }
    LOG(DEBUG) << "parsing of number [" << numberFileContents << "] from [" << s.substr(0,20) << "] failed: " << e.what();
  }
  return number;
}

void removeVaryingContent(std::string &contents)
{
  // in paraview output files there is always a comment "<!-- ... -->" with varying content like timestamp. It is different each
  // run, therefore you cannot compare it to reference values. This functions removes this contents.
  while(contents.find("<!--") != std::string::npos)
  {
    std::size_t pos0 = contents.find("<!--");
    std::size_t pos1 = contents.find("-->", pos0);
    contents = contents.substr(0, pos0) + contents.substr(pos1+3);
  }

  // also remove timestamp and meta information in python files
  // example: {"version": "opendihu 0.1, build Jan 18 2019 20:51:36, C++ 201402, GCC 7.3.0", "meta": "current time: 2019/01/18 20:57:17, hostname: lapsgs05, n ranks: 4", "meshType" ...

  while(contents.find("\"version") != std::string::npos)
  {
    std::size_t pos = contents.find("\"version") + 1;
    std::size_t startPos = pos;
    for (int i = 0; i < 7; i++)
    {
      pos = contents.find("\"", pos)+1;
    }
    pos = contents.find(",", pos)+1;
    pos = contents.find("\"", pos);
    contents = contents.substr(0, startPos-1) + contents.substr(pos);
  }
  while(contents.find("'version") != std::string::npos)
  {
    std::size_t pos = contents.find("'version") + 1;
    std::size_t startPos = pos;
    for (int i = 0; i < 7; i++)
    {
      pos = contents.find("'", pos)+1;
    }
    pos = contents.find(",", pos)+1;
    pos = contents.find("'", pos);
    contents = contents.substr(0, startPos-1) + contents.substr(pos);
  }
}

void assertFileMatchesContent(std::string filename, std::string referenceContents, std::string referenceContents2)
{
  // read in generated exnode file 
  std::ifstream file;

  file.open(filename);
  for (int i = 0; !file.is_open() && i < 5; i++)
  {
    // wait a bit for the files being ready to be opened
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    file.open(filename);
  }

  LOG(DEBUG) << "assertFileMatchesContent: Parse file \"" << filename << "\"";
  ASSERT_TRUE(file.is_open()) << "could not open output file"; 
  
  // reserve memory of size of file
  file.seekg(0, std::ios::end);   
  size_t fileSize = file.tellg();
  std::string fileContents(fileSize, ' ');
  
  // reset file pointer
  file.seekg(0, std::ios::beg);
  
  // read in file contents
  file.read(&fileContents[0], fileSize);

  // remove parts of the contents that are dependend on version or time
  removeVaryingContent(fileContents);
  removeVaryingContent(referenceContents);
  removeVaryingContent(referenceContents2);

  bool referenceContentMatches = true;
  bool referenceContent2Matches = true;
  std::stringstream msg;
  
  // iterate through file contents and compare against referenceContents
  std::string::iterator iterReferenceContents = referenceContents.begin(); 
  for (std::string::iterator iterFileContents = fileContents.begin(); iterFileContents != fileContents.end() && iterReferenceContents != referenceContents.end();)
  {
    //VLOG(1) << "[" << *iterFileContents << "] ?= [" << *iterReferenceContents << "]";
    if (isdigit(*iterFileContents) || *iterFileContents == '.' || *iterFileContents == '-'
      || isdigit(*iterReferenceContents) || *iterReferenceContents == '.' || *iterReferenceContents == '-'
    )   // characters with which a number can start (e.g. ".5", "-1")
    {
      double numberFileContents = parseNumber(iterFileContents, fileContents.end());
      double numberReferenceContents = parseNumber(iterReferenceContents, referenceContents.end());
      if (fabs(numberFileContents - numberReferenceContents) >= 1e-13)
      {
        referenceContentMatches = false;
        msg << "file content of file \"" << filename << "\" is different (referenceContents). The number " << numberFileContents << " in the file does not match reference value " << numberReferenceContents << ", difference " << fabs(numberFileContents - numberReferenceContents) << std::endl;
      }
    }
    else
    {
      // Do not fail if one of the files has more whitespace at a certain position, this can happen if there are different numbers like " -1e-18" and "  1e-18" with different number of whitespaces.
      if (*iterFileContents == ' ' && *iterReferenceContents != ' ')
      {
        iterFileContents++;
        continue;
      }
      if (*iterFileContents != ' ' && *iterReferenceContents == ' ')
      {
        iterFileContents++;
        continue;
      }

      // if character is diferrent
      if(*iterFileContents != *iterReferenceContents)
      {
        msg << "mismatch at character, file: [" << *iterFileContents << "] != reference: [" << *iterReferenceContents << "], pos: " << std::distance(fileContents.begin(),iterFileContents) << " ";
        referenceContentMatches = false;
        //VLOG(1) << "mismatch!";
      }
      iterFileContents++;
      iterReferenceContents++;
    }
  }

  //VLOG(1) << "referenceContentMatches: " << referenceContentMatches;
  
  if (referenceContents2 != "-")
  {
    // iterate through file contents and compare against referenceContents2
    iterReferenceContents = referenceContents2.begin(); 
    for (std::string::iterator iterFileContents = fileContents.begin(); iterFileContents != fileContents.end() && iterReferenceContents != referenceContents2.end();)
    {
      if (isdigit(*iterFileContents) || *iterFileContents == '.' || *iterFileContents == '-'
        || isdigit(*iterReferenceContents) || *iterReferenceContents == '.' || *iterReferenceContents == '-')
      {
        double numberFileContents = parseNumber(iterFileContents, fileContents.end());
        double numberReferenceContents = parseNumber(iterReferenceContents, referenceContents2.end());
        if (fabs(numberFileContents - numberReferenceContents) >= 1e-13) 
        {
          referenceContent2Matches = false;
          msg << "file content of file \"" << filename << "\" is different (referenceContents2). The number " << numberFileContents << " in the file does not match reference value " << numberReferenceContents << ", difference " << fabs(numberFileContents - numberReferenceContents) << std::endl;
        }
      }
      else
      {
        // Do not fail if one of the files has more whitespace at a certain position, this can happen if there are different numbers like " -1e-18" and "  1e-18" with different number of whitespaces.
        if (*iterFileContents == ' ' && *iterReferenceContents != ' ')
        {
          iterFileContents++;
          continue;
        }
        if (*iterFileContents != ' ' && *iterReferenceContents == ' ')
        {
          iterFileContents++;
          continue;
        }

        // if character is diferrent
        if (*iterFileContents, *iterReferenceContents)
        {
          msg << "mismatch at character file: [" << *iterFileContents << "] != reference: [" << *iterReferenceContents << "], pos: " << std::distance(fileContents.begin(),iterFileContents);
          referenceContent2Matches = false;
        }
        iterFileContents++;
        iterReferenceContents++;
      }
    }
  }
  
  //VLOG(1) << "referenceContent2Matches: " << referenceContent2Matches;

  if (!referenceContentMatches && (referenceContents2 == "-" || !referenceContent2Matches))
  {
    if (!referenceContentMatches)
      LOG(INFO) << "file content of file \"" << filename << "\" is different (referenceContents). fileContents: " << std::endl << fileContents << std::endl << ", referenceContents: " << std::endl << referenceContents;
    if (!referenceContent2Matches)
      LOG(INFO) << "file content of file \"" << filename << "\" is different (referenceContents2). fileContents: " << std::endl << fileContents << std::endl << ", referenceContents2: " << std::endl << referenceContents2;
    
    LOG(INFO) << msg.str();
    ASSERT_TRUE(false) << "neither referenceContent nor referenceContent2 matches! " << msg.str();
  }
  
}

void assertParallelEqualsSerialOutputFiles(std::vector<std::string> &outputFilesToCheck)
{
  std::this_thread::sleep_for(std::chrono::milliseconds(100));  // pause execution, such that output files can be closed

  // command line version
  std::stringstream command;
  command << "../../../dependencies/python/install/bin/python3 ../../../scripts/validate_parallel.py ";
  for (auto filename : outputFilesToCheck)
  {
    command << filename << " ";
  }
  int returnValue = system(command.str().c_str());
  ASSERT_EQ(returnValue, 0);
  return;

  // read in python code
  std::ifstream validationScriptFile("../../../scripts/validate_parallel.py");
  ASSERT_TRUE(validationScriptFile.is_open()) << "Could not load validation script";

  std::stringstream source;
  source << validationScriptFile.rdbuf();
  validationScriptFile.close();

  // prepare command line arguments
  int nCommandLineArguments = outputFilesToCheck.size() + 1;
  std::vector<wchar_t *> argvWChar(nCommandLineArguments);
  argvWChar[0] = Py_DecodeLocale(std::string("validate_parallel.py").c_str(), NULL);
  for (int i = 0; i < outputFilesToCheck.size(); i++)
  {
    argvWChar[i+1] = Py_DecodeLocale(outputFilesToCheck[i].c_str(), NULL);
  }

  PySys_SetArgvEx(nCommandLineArguments, argvWChar.data(), 0);

  // wait a bit until files are ready to be opened by the python script
  std::this_thread::sleep_for(std::chrono::milliseconds(20));

  // execute python script
  int ret = PyRun_SimpleString(source.str().c_str());
  if (ret != 0)
    LOG(INFO) << "return value of PyRun_SimpleString is " << ret;

  if (PyErr_Occurred())
  {
    PyErr_Print();
  }

  // load main module and extract config
  PyObject *mainModule = PyImport_AddModule("__main__");
  PyObject *resultVariable = PyObject_GetAttrString(mainModule, "all_tests_successful");

  ASSERT_EQ (resultVariable, Py_True) << "Parallel and serial output do not match!";
}
