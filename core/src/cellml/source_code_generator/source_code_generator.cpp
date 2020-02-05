#include "cellml/source_code_generator/source_code_generator.h"

#include <Python.h>  // has to be the first included header

#include <vector>
#include <iostream>
#include <sstream>
#include "easylogging++.h"
#include "utility/vector_operators.h"
#include "control/dihu_context.h"
#include "utility/string_utility.h"

#ifdef HAVE_OPENCOR
#include "opencor.h"
#endif

CellMLSourceCodeGenerator::CellMLSourceCodeGenerator() : sourceFileSuffix_(".c")
{

}

void CellMLSourceCodeGenerator::initialize(
  std::string inputFilename, int nInstances, int nStates, int nIntermediates,
  const std::vector<int> &parametersUsedAsIntermediate, const std::vector<int> &parametersUsedAsConstant,
  const std::vector<double> &parametersInitialValues
)
{
  sourceFilename_ = inputFilename;
  nInstances_ = nInstances;
  nStates_ = nStates;
  nIntermediates_ = nIntermediates;

  stateNames_.resize(nStates);
  intermediateNames_.resize(nIntermediates);
  statesInitialValues_.resize(nStates);

  parametersUsedAsIntermediate_.assign(parametersUsedAsIntermediate.begin(), parametersUsedAsIntermediate.end());
  parametersUsedAsConstant_.assign(parametersUsedAsConstant.begin(), parametersUsedAsConstant.end());

  nParameters_ = parametersUsedAsIntermediate_.size() + parametersUsedAsConstant_.size();
  parameters_.resize(nParameters_*nInstances_);

  // convert a xml file to a c file using OpenCOR, if necessary
  this->convertFromXmlToC();

  // get initial values from source file and parse source code
  this->parseSourceCodeFile();   // this sets nIntermediatesInSource_

  // set initial values of parameters
  VLOG(1) << ", parameters_.size(): " << parameters_.size();
  if (parametersInitialValues.size() == parameters_.size())
  {
    std::copy(parametersInitialValues.begin(), parametersInitialValues.end(), parameters_.begin());
    LOG(DEBUG) << "parameters size is matching for all instances";
  }
  else
  {
    LOG(DEBUG) << "copy parameters which were given only for one instance to all instances";
    for (int instanceNo=0; instanceNo<nInstances_; instanceNo++)
    {
      for (int j=0; j<nParameters_; j++)
      {
        parameters_[j*nInstances_ + instanceNo] = parametersInitialValues[j];
      }
    }
  }
  VLOG(1) << "parameters_: " << parameters_;
}

void CellMLSourceCodeGenerator::convertFromXmlToC()
{
  if (DihuContext::ownRankNoCommWorld() == 0)
  {
    // open input file
    std::ifstream file(sourceFilename_);

    if (file.is_open())
    {
      std::string line;
      std::getline(file, line);

      if (line.find("<?xml") != std::string::npos)
      {
        // file is an xml file

        // create a c filename
        std::stringstream s;
        s << "src/" << StringUtility::extractBasename(sourceFilename_) << ".c";

        // create src directory
        int ret = system("mkdir -p src");
        if (ret != 0)
        {
          LOG(ERROR) << "Could not create \"src\" directory.";
        }

        std::string cFilename = s.str();

        // call OpenCOR, if available
#ifdef HAVE_OPENCOR
        std::stringstream command;
        command << "\"" << OPENCOR_BINARY << "\" -c CellMLTools::export \"" << sourceFilename_ << "\" \"" << OPENCOR_FORMATS_DIRECTORY << "/C.xml\" > \"" << cFilename << "\"";
        LOG(DEBUG) << "use opencor to convert file: " << command.str();
        ret = system(command.str().c_str());
        if (ret != 0)
        {
          LOG(ERROR) << "Could not convert CellML XML file to c file using OpenCOR: " << command.str();
        }
        else
        {
          sourceFilename_ = cFilename;
        }
#else
        LOG(FATAL) << "OpenCOR is not available to convert XML to c files, but CellMLAdapter has an XML file as input: \"" << sourceFilename_ << "\".";
#endif
      }
    }
  }
}

std::vector<double> &CellMLSourceCodeGenerator::statesInitialValues()
{
  return statesInitialValues_;
}

const std::vector<std::string> &CellMLSourceCodeGenerator::intermediateNames() const
{
  return intermediateNames_;
}

const std::vector<std::string> &CellMLSourceCodeGenerator::stateNames() const
{
  return stateNames_;
}

const int CellMLSourceCodeGenerator::nParameters() const
{
  return nParameters_;
}

std::vector<double> &CellMLSourceCodeGenerator::parameters()
{
  return parameters_;
}

const std::string CellMLSourceCodeGenerator::sourceFilename() const
{
  return sourceFilename_;
}

std::string CellMLSourceCodeGenerator::additionalCompileFlags() const
{
  return additionalCompileFlags_;
}

std::string CellMLSourceCodeGenerator::sourceFileSuffix() const
{
  return sourceFileSuffix_;
}

std::string CellMLSourceCodeGenerator::compilerCommand() const
{
  return compilerCommand_;
}

void CellMLSourceCodeGenerator::generateSourceFile(std::string outputFilename, std::string optimizationType, bool approximateExponentialFunction)
{
  if (optimizationType == "vc")
  {
    generateSourceFileVc(outputFilename, approximateExponentialFunction);
  }
  else if (optimizationType == "simd")
  {
    generateSourceFileSimd(outputFilename);
  }
  else if (optimizationType == "openmp")
  {
    generateSourceFileOpenMP(outputFilename);
  }
}
