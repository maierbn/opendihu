#include "cellml/source_code_generator/00_source_code_generator_base.h"

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

CellmlSourceCodeGeneratorBase::CellmlSourceCodeGeneratorBase() :
  sourceFileSuffix_(".c")
{

}

void CellmlSourceCodeGeneratorBase::initializeNames(std::string inputFilename, int nInstances, int nStates, int nIntermediates)
{
  sourceFilename_ = inputFilename;
  nInstances_ = nInstances;
  nStates_ = nStates;
  nIntermediates_ = nIntermediates;

  stateNames_.resize(nStates);
  intermediateNames_.resize(nIntermediates);
  statesInitialValues_.resize(nStates);

  // convert a xml file to a c file using OpenCOR, if necessary
  this->convertFromXmlToC();

  // get initial values from source file and parse source code
  this->parseNamesInSourceCodeFile();   // this sets nIntermediatesInSource_
}

void CellmlSourceCodeGeneratorBase::initializeSourceCode(
  const std::vector<int> &parametersUsedAsIntermediate, const std::vector<int> &parametersUsedAsConstant,
  std::vector<double> &parametersInitialValues, int maximumNumberOfParameters, double *parameterValues
)
{
  // parametersInitialValues is the list of initial parameter values as given in the settings
  
  parametersUsedAsIntermediate_.assign(parametersUsedAsIntermediate.begin(), parametersUsedAsIntermediate.end());
  parametersUsedAsConstant_.assign(parametersUsedAsConstant.begin(), parametersUsedAsConstant.end());

  nParameters_ = parametersUsedAsIntermediate_.size() + parametersUsedAsConstant_.size();

  // nParameters_ is the number of parameters per instance as determined from mappings

  if (nParameters_ > maximumNumberOfParameters)
  {
    LOG(FATAL) << "There can only be as many parameters as there are intermediates. This is an arbitrary restriction, if you need more parameters, try increasing the number of intermediates in the C++ source file."
      << "Now you have (" << nParameters_ << " parameters and the maximum possible number is " << maximumNumberOfParameters << ".";
  }

  LOG(DEBUG) << "nParameters_: " << nParameters_ << ", parametersInitialValues.size(): " << parametersInitialValues.size() << ", nIntermediates_: " << nIntermediates_ << ", nInstances_: " << nInstances_;

  // set initial values of parameters
  if (parametersInitialValues.size() == nParameters_)
  {
    LOG(DEBUG) << "parameters size is matching for one instances";
    LOG(DEBUG) << "copy parameters which were given only for one instance to all instances";
    for (int instanceNo = 0; instanceNo < nInstances_; instanceNo++)
    {
      for (int j = 0; j < nParameters_; j++)
      {
        // parameterValues has struct of array memory layout with space for a total of nIntermediates_ parameters [i0p0, i1p0, i2p0, ... i0p1, i1p1, i2p1, ...]
        parameterValues[j*nInstances_ + instanceNo] = parametersInitialValues[j];
        
        VLOG(1) << "  " << instanceNo << "," << j << " set index " << (j*nInstances_ + instanceNo) << " to value " << parametersInitialValues[j] << "=" << parameterValues[j*nInstances_ + instanceNo];
      }
    }
  }
  else if (parametersInitialValues.size() == nParameters_ * nInstances_)
  {
    LOG(DEBUG) << "parameters size is matching for all instances";
    for (int instanceNo = 0; instanceNo < nInstances_; instanceNo++)
    {
      for (int j = 0; j < nParameters_; j++)
      {
        // parameters are given in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
        parameterValues[j*nInstances_ + instanceNo] = parametersInitialValues[instanceNo*nParameters_ + j];
      }
    }
  }
  else 
  {
    LOG(WARNING) << "In CellML: There should be " << nParameters_ << " parameters but " << parametersInitialValues.size()
      << " initial values are given by \"parametersInitialValues\". Using default values 0.";
    parametersInitialValues.resize(nParameters_, 0.0);
  
    LOG(DEBUG) << "copy parameters which were given only for one instance to all instances";
    for (int instanceNo = 0; instanceNo < nInstances_; instanceNo++)
    {
      for (int j = 0; j < nParameters_; j++)
      {
        // parameterValues has struct of array memory layout with space for a total of nIntermediates_ parameters [i0p0, i1p0, i2p0, ... i0p1, i1p1, i2p1, ...]
        parameterValues[j*nInstances_ + instanceNo] = parametersInitialValues[j];
      }
    }
  }
  
#ifndef NDEBUG
  std::stringstream s;
  for (int instanceNo = 0; instanceNo < nInstances_; instanceNo++)
  {
    for (int j = 0; j < nParameters_; j++)
    {
      // parameters are given in array of struct ordering: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
      s << parameterValues[j*nInstances_ + instanceNo] << " ";
    }
    s << ",";
  }
  LOG(DEBUG) << "parameterValues (only up to nParameters=" << nParameters_ << " parameter, allocated space is for " << nIntermediates_ << " parameters: \n" << s.str();
#endif  
  
  // parse all the source code from the model file
  this->parseSourceCodeFile();

  // Generate the rhs code for a single instance. This is needed for computing the equilibrium of the states.
  this->generateSingleInstanceCode();
}

void CellmlSourceCodeGeneratorBase::convertFromXmlToC()
{
  std::string cFilename = sourceFilename_;

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
      cFilename = s.str();
    }
  }

  // check if conversion is required, this is performed on rank 0
  bool useOpenCor = false;
  if (DihuContext::ownRankNoCommWorld() == 0)
  {
    // check if file already exists
    std::ifstream cFile(cFilename.c_str());
    if (cFile.is_open())
    {
      LOG(DEBUG) << "C file \"" << cFilename << "\" already exists.";
      sourceFilename_ = cFilename;
    }
    else
    {

      // create src directory
      int ret = system("mkdir -p src");
      if (ret != 0)
      {
        LOG(ERROR) << "Could not create \"src\" directory.";
      }

          // should call OpenCOR later, if available
#ifdef HAVE_OPENCOR
      useOpenCor = true;
#else
      LOG(FATAL) << "OpenCOR is not available to convert XML to c files, but CellMLAdapter has an XML file as input: \"" << sourceFilename_ << "\".";
#endif
    }
  }

  LOG(DEBUG) << "useOpenCor=" << useOpenCor << ", cFilename: " << cFilename << ", sourceFilename_: " << sourceFilename_;

  if (useOpenCor)
  {
    // do conversion on rank 0
    if (DihuContext::ownRankNoCommWorld() == 0)
    {
      std::stringstream command;
      command << "\"" << OPENCOR_BINARY << "\" -c CellMLTools::export \"" << sourceFilename_ << "\" \"" << OPENCOR_FORMATS_DIRECTORY << "/C.xml\" > \"" << cFilename << "\"";
      LOG(DEBUG) << "Use opencor to convert file: \n" << command.str();
      int ret = system(command.str().c_str());
      if (ret != 0)
      {
        LOG(FATAL) << "Could not convert CellML XML file \"" << sourceFilename_ << "\" to c file \"" << cFilename << "\" using OpenCOR: \n" << command.str();
      }
    }

    // wait on all ranks until conversion is finished
    MPIUtility::handleReturnValue(MPI_Barrier(DihuContext::partitionManager()->rankSubsetForCollectiveOperations()->mpiCommunicator()), "MPI_Barrier");
  }

  sourceFilename_ = cFilename;
}

void CellmlSourceCodeGeneratorBase::generateSingleInstanceCode()
{

  LOG(DEBUG) << "generateSingleInstanceCode";

  std::stringstream sourceCode;

  sourceCode << std::endl
    << "// compute the rhs for a single instance, this can be used for computation of the equilibrium values of the states" << std::endl
    << "#ifdef __cplusplus\n"
    << "extern \"C\"\n"
    << "#endif\n"
    << "void computeCellMLRightHandSideSingleInstance("
    << "void *context, double t, double *states, double *rates, double *intermediates, double *parameters)" << std::endl << "{" << std::endl;

  sourceCode << "  double VOI = t;   /* current simulation time */" << std::endl;
  sourceCode << std::endl << "  /* define constants */" << std::endl
    << "  double CONSTANTS[" << this->nConstants_ << "];" << std::endl;

  // add assignments of constant values
  for (std::string constantAssignmentsLine : constantAssignments_)
  {
    sourceCode << "  " << constantAssignmentsLine << std::endl;
  }
  sourceCode << std::endl;

  // loop over lines of cellml code
  for (code_expression_t &codeExpression : cellMLCode_.lines)
  {

    if (codeExpression.type != code_expression_t::commented_out)
    {
      sourceCode << "  ";
      codeExpression.visitLeafs([&sourceCode,this](code_expression_t &expression, bool isFirstVariable)
      {
        switch(expression.type)
        {
        case code_expression_t::variableName:
          sourceCode << expression.code << "[" << expression.arrayIndex<< "]";
          break;

        case code_expression_t::otherCode:
          sourceCode << expression.code;
          break;

        case code_expression_t::commented_out:
          sourceCode << "  // (not assigning to a parameter) " << expression.code;
          break;

        default:
          break;
        };
      });

      sourceCode << std::endl;
    }
  }

  sourceCode << "}\n";

  singleInstanceCode_ = sourceCode.str();
}

std::vector<double> &CellmlSourceCodeGeneratorBase::statesInitialValues()
{
  return statesInitialValues_;
}

const std::vector<std::string> &CellmlSourceCodeGeneratorBase::intermediateNames() const
{
  return intermediateNames_;
}

const std::vector<std::string> &CellmlSourceCodeGeneratorBase::stateNames() const
{
  return stateNames_;
}

const std::vector<std::string> &CellmlSourceCodeGeneratorBase::constantNames() const
{
  return constantNames_;
}

const int CellmlSourceCodeGeneratorBase::nParameters() const
{
  return nParameters_;
}

const std::string CellmlSourceCodeGeneratorBase::sourceFilename() const
{
  return sourceFilename_;
}

std::string CellmlSourceCodeGeneratorBase::additionalCompileFlags() const
{
  return additionalCompileFlags_;
}

std::string CellmlSourceCodeGeneratorBase::sourceFileSuffix() const
{
  return sourceFileSuffix_;
}

std::string CellmlSourceCodeGeneratorBase::compilerCommand() const
{
  return compilerCommand_;
}
