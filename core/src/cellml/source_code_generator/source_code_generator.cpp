#include "cellml/source_code_generator/source_code_generator.h"

#include <Python.h>  // has to be the first included header

#include <vector>
#include <iostream>
#include "easylogging++.h"
#include "utility/vector_operators.h"

CellMLSourceCodeGenerator::CellMLSourceCodeGenerator()
{

}

void CellMLSourceCodeGenerator::initialize(
  std::string inputFilename, int nInstances, int nStates, int nIntermediates,
  const std::vector<int> &parametersUsedAsIntermediate, const std::vector<int> &parametersUsedAsConstant,
  const std::vector<double> &parametersInitialValues
)
{
  nInstances_ = nInstances;

  stateNames_.resize(nStates);
  intermediateNames_.resize(nIntermediates);
  statesInitialValues_.resize(nStates);

  parametersUsedAsIntermediate_.assign(parametersUsedAsIntermediate.begin(), parametersUsedAsIntermediate.end());
  parametersUsedAsConstant_.assign(parametersUsedAsConstant.begin(), parametersUsedAsConstant.end());

  nParameters_ = parametersUsedAsIntermediate_.size() + parametersUsedAsConstant_.size();
  parameters_.resize(nParameters_*nInstances_);

  // get initial values from source file and parse source code
  this->parseSourceCodeFile();   // this sets nIntermediatesInSource_

  // check number of intermediates in source file
  if (nIntermediatesInSource_ > nIntermediates_)
  {
    LOG(FATAL) << "CellML source file needs " << nIntermediatesInSource_ << " intermediates, but CellMLAdapter only supports " << nIntermediates_
      << ". You have to set the correct number in the c++ file and recompile." << std::endl
      << "(Use \"CellMLAdapter<" << nStates << "," << nIntermediatesInSource_ << ">\".)";
  }
  else if (nIntermediatesInSource_ != nIntermediates_)
  {
    LOG(WARNING) << "CellML source file needs " << nIntermediatesInSource_ << " intermediates, and CellMLAdapter supports " << nIntermediates_
      << ". You should recompile with the correct number to avoid performance penalties.";
  }


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

/*
bool CellMLSourceCodeGenerator::scanSourceFile(std::string sourceFilename)
{
  LOG(TRACE) << "scanSourceFile";

  // parse source file, set initial values for states (only one instance) and nParameters_, nConstants_ and nIntermediatesInSource_
  nIntermediatesInSource_ = 0;

  // read in source from file
  std::ifstream sourceFile(sourceFilename.c_str());
  if (!sourceFile.is_open())
  {
    LOG(WARNING) << "Could not open source file \"" << sourceFilename << "\" for reading initial values.";
    return false;
  }
  else
  {
    std::string name;  // the parsed name of a specifier that follows

    // read whole file contents
    std::stringstream source;
    source << sourceFile.rdbuf();
    sourceFile.close();

    // step through lines of simd file
    while(!source.eof())
    {
      std::string line;
      getline(source, line);

      if (line.find(" * STATES") == 0)  // line in OpenCOR generated input file of type " * STATES[55] is P_C_SR in component razumova (milliM)."
      {
        // parse name of state
        unsigned int index = atoi(line.substr(10,line.find("]")-10).c_str());
        int posBegin = line.find("is",12)+3;
        int posEnd = line.rfind(" in");
        name = line.substr(posBegin,posEnd-posBegin);
        LOG(DEBUG) << "index= " << index << ", this->stateNames_.size() = " << this->stateNames_.size();
        if (index >= this->stateNames_.size())
        {
          LOG(FATAL) << "The CellML file \"" << sourceFilename << "\" contains more than " << index << " states "
            << " but only " << this->stateNames_.size() << " were given as template argument to CellMLAdapter.";
        }
        this->stateNames_[index] = name;
      }
      else if (line.find(" * ALGEBRAIC") == 0)  // line in OpenCOR generated input file of type " * ALGEBRAIC[35] is g_Cl in component sarco_Cl_channel (milliS_per_cm2)."
      {
        // parse name of intermediate
        unsigned int index = atoi(line.substr(13,line.find("]")-13).c_str());
        int posBegin = line.find("is",15)+3;
        int posEnd = line.rfind(" in");
        name = line.substr(posBegin,posEnd-posBegin);
        LOG(DEBUG) << "index= " << index << ", this->intermediateNames_.size() = " << this->intermediateNames_.size();
        if (index >= this->intermediateNames_.size())
        {
          LOG(FATAL) << "The CellML file \"" << sourceFilename << "\" contains more than " << index << " intermediates "
            << " but only " << this->intermediateNames_.size() << " were given as template argument to CellMLAdapter.";
        }
        this->intermediateNames_[index] = name;
      }
      else if (line.find("STATES[") == 0)   // line contains assignment in OpenCOR generated input file
      {
        // parse initial value of state
        unsigned int index = atoi(line.substr(7,line.find("]",7)-7).c_str());
        double value = atof(line.substr(line.find("= ")+2).c_str());
        statesInitialValues_[index] = value;
      }
      else if (line.find("ALGEBRAIC[") == 0)  // assignment to an algebraic variable in both OpenCMISS and OpenCOR generated files, in OpenCMISS generated files, this does not count towards the algebraic variables that are hold by opendihu
      {
        int algebraicIndex = atoi(line.substr(10,line.find("]",10)-10).c_str());
        nIntermediatesInSource_ = std::max(nIntermediatesInSource_, algebraicIndex+1);
      }
      else if (line.find("CONSTANTS[") != std::string::npos)  // usage of a constant
      {
        std::string substr(line.substr(line.find("CONSTANTS[")+10,line.find("]",line.find("CONSTANTS[")+10)-line.find("CONSTANTS[")-10));
        int index = atoi(substr.c_str());
        this->nConstants_ = std::max(this->nConstants_, index+1);
      }

      // if the rhs routine is reached in the OpenCOR file, there are no more initializations, there quit processing of the source file
      if (line.find("computeCellMLRightHandSide") != std::string::npos)
      {
        break;
      }
    }
  }

  return true;
}
*/
void CellMLSourceCodeGenerator::generateSourceFile(std::string outputFilename, std::string optimizationType) const
{
  if (optimizationType == "vc")
  {
    generateSourceFileExplicitVectorization(outputFilename);
  }
  else if (optimizationType == "simd")
  {
    generateSourceFileSimdAutovectorization(outputFilename);
  }
  else if (optimizationType == "openmp")
  {

  }
}
