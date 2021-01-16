#include "cellml/source_code_generator/04_generator_gpu.h"

#include <Python.h>  // has to be the first included header

#include "utility/string_utility.h"
#include "output_writer/generic.h"

#include <vector>
#include <iostream>
#include <Vc/Vc>
#include "easylogging++.h"

void CellmlSourceCodeGeneratorGpu::
generateSourceFileGpu(std::string outputFilename)
{
  std::stringstream sourceCode;
  sourceCode << "#include <math.h>" << std::endl
    << "#include <omp.h>" << std::endl
    << cellMLCode_.header << std::endl;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  sourceCode << std::endl << "/* This function was created by opendihu at " << StringUtility::timeToString(&tm)  //std::put_time(&tm, "%d/%m/%Y %H:%M:%S")
    << ".\n * It is designed for " << this->nInstances_ << " instances of the CellML problem.\n "
    << " * The \"optimizationType\" is \"gpu\". (Other options are \"vc\", \"openmp\" and \"simd\".) */" << std::endl
    << "void computeCellMLRightHandSide("
    << "void *context, double t, double *states, double *rates, double *algebraics, double *parameters)" << std::endl << "{" << std::endl;

  sourceCode << "  double VOI = t;   /* current simulation time */" << std::endl;
  sourceCode << std::endl << "  /* define constants */" << std::endl
    << "  double CONSTANTS[" << this->nConstants_ << "];" << std::endl;

  // add assignments of constant values
  for (std::string constantAssignmentsLine : constantAssignments_)
  {
    sourceCode << "  " << constantAssignmentsLine << std::endl;
  }

  sourceCode << std::endl
    << "  #pragma omp target parallel for map(to:states,t,parameters) map(from:rates,algebraics)" << std::endl
    << "  for (int i = 0; i < " << this->nInstances_ << "; i++)" << std::endl
    << "  {" << std::endl;

  // loop over lines of cellml code
  for (code_expression_t &codeExpression : cellMLCode_.lines)
  {
    if (codeExpression.type != code_expression_t::commented_out)
    {

      sourceCode << "    ";
      codeExpression.visitLeafs([&sourceCode,this](CellmlSourceCodeGeneratorOpenMp::code_expression_t &expression, bool isFirstVariable)
      {

        switch(expression.type)
        {
        case code_expression_t::variableName:

          if (expression.code == "CONSTANTS")
          {
            // constants only exist once for all instances
            sourceCode << expression.code << "[" << expression.arrayIndex<< "]";
          }
          else
          {
            // all other variables (states, rates, algebraics, parameters) exist for every instance
            sourceCode << expression.code << "[" << expression.arrayIndex * this->nInstances_ << "+i]";
          }
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

  VLOG(2) << "write end of for loop (closing })";
  sourceCode << "  }" << std::endl;

  // add footer
  sourceCode << cellMLCode_.footer << std::endl;

  // add code for a single instance
  sourceCode << singleInstanceCode_;

  // write out source file
  std::ofstream sourceCodeFile;
  OutputWriter::Generic::openFile(sourceCodeFile, outputFilename);
  if (!sourceCodeFile.is_open())
  {
    LOG(FATAL) << "Could not write to file \"" << outputFilename << "\".";
  }
  else
  {
    std::string fileContents = sourceCode.str();
    sourceCodeFile << fileContents;
    sourceCodeFile.close();
  }

  additionalCompileFlags_ = "-fopenmp -foffload=\"-O3 -v -lm\"";
  compilerCommand_ = C_COMPILER_COMMAND;
  sourceFileSuffix_ = ".c";
}

void CellmlSourceCodeGeneratorGpu::
generateSourceFastMonodomainGpu(bool approximateExponentialFunction, int nFibersToCompute, int nInstancesToComputePerFiber, 
                                int nParametersPerInstance, bool hasAlgebraicsForTransfer,
                                std::string &headerCode, std::string &mainCode)
{
  std::stringstream sourceCodeHeader;
  std::stringstream sourceCodeMain;
  
  std::set<std::string> helperFunctions;   //< functions found in the CellML code that need to be provided, usually the pow2, pow3, etc. helper functions for pow(..., 2), pow(...,3) etc.

  // replace pow and ?: functions
  const bool useVc = false;
  preprocessCode(helperFunctions, useVc);

  sourceCodeHeader << "#include <cmath>\n"
    << "#include <omp.h>\n"
    << "#include <iostream>\n"
    << "#include <vector>\n\n"
    << cellMLCode_.header << std::endl
    << R"(
real log(real x)
{
  // Taylor expansion of the log function around x=1
  // note: std::log does not work on GPU!
  real t = x-1;
  real t2 = t*t;
  return t - 0.5*t2 + 1./3*t2*t - 0.25*t2*t2;
}
)";
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  sourceCodeHeader << std::endl << "/* This file was created by opendihu at " << StringUtility::timeToString(&tm)  //std::put_time(&tm, "%d/%m/%Y %H:%M:%S")
    << ".\n * It is designed for the FastMonodomainSolver and contains code for offloading to GPU.\n "
    << " */\n";

  VLOG(1) << "call defineHelperFunctions with helperFunctions: " << helperFunctions;

  // define helper functions
  const bool useReal = true;
  sourceCodeHeader << defineHelperFunctions(helperFunctions, approximateExponentialFunction, false, useReal);
    
  // define main code that computes the rhs
  
  long int nInstancesToCompute = nInstancesToComputePerFiber * nFibersToCompute;
  std::string indent(10, ' ');
  
  sourceCodeMain << "\n" << indent << "// CellML define constants\n";
    
  // add assignments of constant values
  for (std::string constantAssignmentsLine : constantAssignments_)
  {
    constantAssignmentsLine = StringUtility::replaceAll(constantAssignmentsLine, "CONSTANTS[", "constant");
    constantAssignmentsLine = StringUtility::replaceAll(constantAssignmentsLine, "]", "");

    sourceCodeMain << indent << "const real " << constantAssignmentsLine << std::endl;
  }
  
  // loop over instances on the current fiber  
  sourceCodeMain << "\n" << indent << "// compute new rates, rhs(y_n)\n";

  // loop over lines of cellml code
  for (code_expression_t &codeExpression : cellMLCode_.lines)
  {
    if (codeExpression.type != code_expression_t::commented_out)
    {
      std::stringstream sourceCodeLine;
      bool isCommentedOut = false;
      bool isFirst = true;
      
      codeExpression.visitLeafs([&sourceCodeLine,&isCommentedOut,nInstancesToCompute,&isFirst,this](
        CellmlSourceCodeGeneratorVc::code_expression_t &expression, bool isFirstVariable)
      {
        switch(expression.type)
        {
          case code_expression_t::variableName:

            if (expression.code == "CONSTANTS")
            {
              // constants only exist once for all instances
              if (isFirst)
                sourceCodeLine << "const real ";
              sourceCodeLine << "constant" << expression.arrayIndex;
            }
            else
            {
              // all other variables (states, rates, algebraics, parameters) exist for every instance
              if (expression.code == "states")
              {
                if (expression.arrayIndex == 0)
                {
                  sourceCodeLine << "vmValues[instanceToComputeNo]";
                }
                else
                {
                  sourceCodeLine << "states[" << expression.arrayIndex*nInstancesToCompute << "+instanceToComputeNo]";
                }
              }
              else if (expression.code == "rates")
              {
                if (isFirst)
                  sourceCodeLine << "const real ";
                sourceCodeLine << "rate" << expression.arrayIndex;
              }
              else if (expression.code == "algebraics")
              {
                if (isFirst)
                  sourceCodeLine << "const real ";
                sourceCodeLine << "algebraic" << expression.arrayIndex;
              }
              else if (expression.code == "parameters")
              {
                sourceCodeLine << "parameters[" << expression.arrayIndex*nInstancesToCompute << "+instanceToComputeNo]";
              }
              else
              {
                LOG(FATAL) << "unhandled variable type \"" << expression.code << "\".";
              }
            }
            break;

          case code_expression_t::otherCode:
            sourceCodeLine << expression.code;
            break;

          case code_expression_t::commented_out:
            sourceCodeLine << "// (not assigning to a parameter) " << expression.code;
            isCommentedOut = true;
            break;

          default:
            break;
        }
        isFirst = false;
      });
      
      sourceCodeMain << indent << sourceCodeLine.str() << std::endl;
    }
  }
  sourceCodeMain << "\n"
    << "          // algebraic step\n"
    << "          // compute y* = y_n + dt*rhs(y_n), y_n = state, rhs(y_n) = rate, y* = intermediateState\n";

  sourceCodeMain << indent << "states[0+instanceToComputeNo] = vmValues[instanceToComputeNo] + dt0D*rate0;\n";

  for (int stateNo = 1; stateNo < this->nStates_; stateNo++)
  {
    sourceCodeMain << indent
       << "states[" << stateNo*nInstancesToCompute << "+instanceToComputeNo] = states[" << stateNo*nInstancesToCompute << "+instanceToComputeNo] + dt0D*rate" << stateNo << ";\n";
  }
  
  sourceCodeMain << "\n"
      << R"(
          // if stimulation, set value of Vm (state0)
          if (stimulateCurrentPoint)
          {
            states[0+instanceToComputeNo] = valueForStimulatedPoint;
          })";
  
  sourceCodeMain << R"(
          // compute new rates, rhs(y*)
)";

  // loop over lines of cellml code
  for (code_expression_t &codeExpression : cellMLCode_.lines)
  {
    if (codeExpression.type != code_expression_t::commented_out)
    {
      std::stringstream sourceCodeLine;
      bool isCommentedOut = false;
      bool isFirst = true;
      
      codeExpression.visitLeafs([&sourceCodeLine,&isCommentedOut,nInstancesToCompute,&isFirst,this](
        CellmlSourceCodeGeneratorVc::code_expression_t &expression, bool isFirstVariable)
      {
        switch(expression.type)
        {
          case code_expression_t::variableName:

            if (expression.code == "CONSTANTS")
            {
              // constants only exist once for all instances
              if (isFirst)
                sourceCodeLine << "const real ";
              sourceCodeLine << "constant" << expression.arrayIndex;
            }
            else
            {
              // all other variables (states, rates, algebraics, parameters) exist for every instance
              if (expression.code == "states")
              {
                sourceCodeLine << "states[" << expression.arrayIndex*nInstancesToCompute << "+instanceToComputeNo]";
              }
              else if (expression.code == "rates")
              {
                if (isFirst)
                  sourceCodeLine << "const real ";
                sourceCodeLine << "intermediateRate" << expression.arrayIndex;
              }
              else if (expression.code == "algebraics")
              {
                if (isFirst)
                  sourceCodeLine << "const real ";
                sourceCodeLine << "intermediateAlgebraic" << expression.arrayIndex;
              }
              else if (expression.code == "parameters")
              {
                sourceCodeLine << "parameters[" << expression.arrayIndex*nInstancesToCompute << "+instanceToComputeNo]";
              }
              else
              {
                LOG(FATAL) << "unhandled variable type \"" << expression.code << "\".";
              }
            }
            break;

          case code_expression_t::otherCode:
            sourceCodeLine << expression.code;
            break;

          case code_expression_t::commented_out:
            sourceCodeLine << "    // (not assigning to a parameter) " << expression.code;
            isCommentedOut = true;
            break;

          default:
            break;
        }
        isFirst = false;
      });
      
      sourceCodeMain << indent << sourceCodeLine.str() << std::endl;
    }
  }

  sourceCodeMain << R"(
          // final step
          // y_n+1 = y_n + 0.5*[rhs(y_n) + rhs(y*)]
)";

  sourceCodeMain << indent << "vmValues[instanceToComputeNo] += 0.5*dt0D*(rate0 + intermediateRate0);\n";
  
  for (int stateNo = 1; stateNo < this->nStates_; stateNo++)
  {
    sourceCodeMain << indent << "states[" << stateNo*nInstancesToCompute << "+instanceToComputeNo] += 0.5*dt0D*(rate" << stateNo << " + intermediateRate" << stateNo << ");\n";
  }

  sourceCodeMain << R"(
          if (stimulateCurrentPoint)
          {
            vmValues[instanceToComputeNo] = valueForStimulatedPoint;
          }  
  )";
  sourceCodeMain << R"(
          // store algebraics for transfer
          if (storeAlgebraicsForTransfer)
          {)";
  if (hasAlgebraicsForTransfer)
  {
    sourceCodeMain << R"(
            for (int i = 0; i < nAlgebraicsForTransferIndices; i++)
            {
              const int algebraicIndex = algebraicsForTransferIndices[i];

              switch (algebraicIndex)
              {
)";

    // loop over algebraics and generate code to copy the updated algebraic values to the algebraics
    for (int algebraicNo = 0; algebraicNo < this->nAlgebraics_; algebraicNo++)
    {
      // only of the algebraic was computed and not replaced by a parameter
      if (std::find(this->parametersUsedAsAlgebraic_.begin(), this->parametersUsedAsAlgebraic_.end(), algebraicNo)
         != this->parametersUsedAsAlgebraic_.end())
      {
        sourceCodeMain << indent << "      // case " << algebraicNo << ": is a parameter\n";
      }
      else
      {
        sourceCodeMain << indent << "      case " << algebraicNo << ":\n"
          << indent << "        algebraicsForTransfer[i*nInstancesToCompute + instanceToComputeNo] = intermediateAlgebraic" << algebraicNo << ";\n"
          << indent << "        break;\n";
      }
    }
    sourceCodeMain << R"(
              }
            })";
  }
  sourceCodeMain << R"(
            
            for (int i = 0; i < nStatesForTransferIndices; i++)
            {
              const int stateIndex = statesForTransferIndices[i];

              switch (stateIndex)
              {
                case 0:
                  statesForTransfer[i*nInstancesToCompute + instanceToComputeNo] = vmValues[instanceToComputeNo];
                  break;
)";

  // loop over states and generate code to copy the updated state values to the statesForTransfer
  for (int stateNo = 1; stateNo < this->nStates_; stateNo++)
  {
    // only of the algebraic was computed and not replaced by a parameter
    sourceCodeMain << indent << "      case " << stateNo << ":\n"
      << indent << "        statesForTransfer[i*nInstancesToCompute + instanceToComputeNo] = states[" << stateNo*nInstancesToCompute << "+instanceToComputeNo];\n"
      << indent << "        break;\n";
  }
  sourceCodeMain << R"(
              }
            }              
          }
)";

  // assign output to reference arguments
  headerCode = sourceCodeHeader.str();
  mainCode = sourceCodeMain.str();
  
  additionalCompileFlags_ = "-fopenmp -foffload=\"-O3 -v -lm\"";
  compilerCommand_ = CXX_COMPILER_COMMAND;
  sourceFileSuffix_ = ".cpp";
}
