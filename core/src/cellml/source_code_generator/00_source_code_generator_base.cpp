#include "cellml/source_code_generator/00_source_code_generator_base.h"

#include <Python.h>  // has to be the first included header

#include <vector>
#include <iostream>
#include <thread>
#include <chrono>
#include <sys/stat.h> // stat
#include <unistd.h>   // stat
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

void CellmlSourceCodeGeneratorBase::initializeNames(std::string inputFilename, int nInstances, int nStates, int nAlgebraics)
{
  sourceFilename_ = inputFilename;
  nInstances_ = nInstances;
  nStates_ = nStates;
  nAlgebraics_ = nAlgebraics;

  stateNames_.resize(nStates);
  algebraicNames_.resize(nAlgebraics);
  statesInitialValues_.resize(nStates);

  // convert a xml file to a c file using OpenCOR, if necessary
  this->convertFromXmlToC();

  // get initial values from source file and parse source code
  this->parseNamesInSourceCodeFile();   // this sets nAlgebraicsInSource_
}

void CellmlSourceCodeGeneratorBase::initializeSourceCode(
  const std::vector<int> &parametersUsedAsAlgebraic, const std::vector<int> &parametersUsedAsConstant,
  std::vector<double> &parametersInitialValues, int maximumNumberOfParameters, double *parameterValues
)
{
  // parametersInitialValues is the list of initial parameter values as given in the settings
  
  parametersUsedAsAlgebraic_.assign(parametersUsedAsAlgebraic.begin(), parametersUsedAsAlgebraic.end());
  parametersUsedAsConstant_.assign(parametersUsedAsConstant.begin(), parametersUsedAsConstant.end());

  nParameters_ = parametersUsedAsAlgebraic_.size() + parametersUsedAsConstant_.size();

  // nParameters_ is the number of parameters per instance as determined from mappings

  if (nParameters_ > maximumNumberOfParameters)
  {
    LOG(FATAL) << "There can only be as many parameters as there are algebraics. This is an arbitrary restriction, if you need more parameters, try increasing the number of algebraics in the C++ source file."
      << "Now you have (" << nParameters_ << " parameters and the maximum possible number is " << maximumNumberOfParameters << ").";
  }

  LOG(DEBUG) << "nParameters_: " << nParameters_ << ", parametersInitialValues.size(): " << parametersInitialValues.size() << ", nAlgebraics_: " << nAlgebraics_ << ", nInstances_: " << nInstances_;

  if (parametersInitialValues.size() != nParameters_ && parametersInitialValues.size() != nParameters_ * nInstances_)
  {
    // the given number of parameters does not match the number of initial values in the settings
    LOG(WARNING) << "In CellML: There should be " << nParameters_
      << " = " << parametersUsedAsAlgebraic_.size() << "+" << parametersUsedAsConstant_.size()
      << " parameters (algebraics: " << parametersUsedAsAlgebraic_ << ", constants: " << parametersUsedAsConstant_ << ") or "
      << nParameters_ * nInstances_ << " parameters in array of struct format, but "  << parametersInitialValues.size()
      << " initial values are given by \"parametersInitialValues\".";
    if (parametersInitialValues.size() < nParameters_)
      LOG(WARNING) << "  Using '0' as default value for " << (nParameters_ - parametersInitialValues.size()) << "undefiend parameters.";
    else
      LOG(WARNING) << "  Using only the first " << nParameters_ << " parameters";

    parametersInitialValues.resize(nParameters_, 0.0);
  }

  // set initial values of parameters
  if (parametersInitialValues.size() == nParameters_)
  {
    // only parameters for one instance: use these for all instances
    LOG(DEBUG) << "parameters size is matching for one instances";
    LOG(DEBUG) << "copy parameters which were given only for one instance to all " << nInstances_ << " instances";
    for (int instanceNo = 0; instanceNo < nInstances_; instanceNo++)
    {
      for (int j = 0; j < nParameters_; j++)
      {
        // parameterValues has struct of array memory layout with space for a total of nAlgebraics_ parameters [i0p0, i1p0, i2p0, ... i0p1, i1p1, i2p1, ...]
        parameterValues[j*nInstances_ + instanceNo] = parametersInitialValues[j];
        
        VLOG(1) << "  " << instanceNo << "," << j << " set index " << (j*nInstances_ + instanceNo) << " to value " << parametersInitialValues[j] << "=" << parameterValues[j*nInstances_ + instanceNo];
      }
    }
  }
  else // parametersInitialValues.size() == nParameters_ * nInstances_
  {
    // parameters are given for all instances
    // contiguous array of struct: inst1p0, inst1p1, ... instNp0, instNp1... instNpM
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

  // assemble information about parameter mappings
  std::stringstream s;
  bool firstEntry = true;
  int parameterIndex = 0;
  for (int algebraicIndex : parametersUsedAsAlgebraic_)
  {
    if (!firstEntry)
      s << ", ";
    s << "\n  parameter " << parameterIndex << " maps to \"" << algebraicNames_[algebraicIndex] << "\" ("
      << "ALGEBRAIC[" << algebraicIndex << "]), "
      << "initial value: " << parameterValues[parameterIndex*nInstances_] << "";

    firstEntry = false;
    parameterIndex++;
  }

  for (int constantIndex : parametersUsedAsConstant_)
  {
    if (!firstEntry)
      s << ", ";
    s << "\n  parameter " << parameterIndex << " maps to \"" << constantNames_[constantIndex] << "\" ("
      << "CONSTANTS[" << constantIndex << "]), "
      << "initial value: " << parameterValues[parameterIndex*nInstances_] << "";

    firstEntry = false;
    parameterIndex++;
  }

  std::stringstream cellmlMessage;
  cellmlMessage << "CellML file \"" << sourceFilename_ << "\" with "
    << nStates_ << " state" << (nStates_!=1? "s" : "") <<", " << nAlgebraicsInSource_
    << " algebraic" << (nAlgebraicsInSource_!=1? "s" : "")
    <<", specified " << nParameters_ << " parameter" << (nParameters_!=1? "s" : "") << ": " << s.str() << "\n";

  // only print message if it has not already been printed
  static std::vector<std::string> cellmlMessages;
  if (std::find(cellmlMessages.begin(), cellmlMessages.end(), cellmlMessage.str()) == cellmlMessages.end())
  {
    LOG(INFO) << cellmlMessage.str();
    cellmlMessages.push_back(cellmlMessage.str());
  }

#ifndef NDEBUG
  std::stringstream message;
  for (int instanceNo = 0; instanceNo < nInstances_; instanceNo++)
  {
    for (int j = 0; j < nParameters_; j++)
    {
      // parameters are given in array of struct ordering: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
      message << parameterValues[j*nInstances_ + instanceNo] << " ";
    }
    message << ",";
  }
  LOG(DEBUG) << "parameterValues (only up to nParameters=" << nParameters_ << " parameter, allocated space is for " << nAlgebraics_ << " parameters: \n" << message.str();
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
  else
  {
    LOG(FATAL) << "Could not open CellML file \"" << sourceFilename_ << "\".";
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
  }

  //LOG(DEBUG) << "MPI barrier, wait on all ranks until conversion is finished, n ranks: " << DihuContext::partitionManager()->rankSubsetForCollectiveOperations()->size();

  // wait on all ranks until conversion is finished, this fails with multiple ranks in x,y and z direction
  //MPIUtility::handleReturnValue(MPI_Barrier(DihuContext::partitionManager()->rankSubsetForCollectiveOperations()->mpiCommunicator()), "MPI_Barrier");

  // wait until file exists
  for (;;)
  {
    struct stat buffer;
    // if file cFilename exists, exit infinite loop
    if (stat(cFilename.c_str(), &buffer) == 0)
    {
      break;
    }

    // yield to other processes
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    std::this_thread::yield();
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
    << "void *context, double t, double *states, double *rates, double *algebraics, double *parameters)" << std::endl << "{" << std::endl;

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

const std::vector<std::string> &CellmlSourceCodeGeneratorBase::algebraicNames() const
{
  return algebraicNames_;
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
