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
    << "#include <vector>\n"
    << "#pragma omp declare target\n\n"
    << cellMLCode_.header << std::endl
    << R"(
double log(double x)
{
  return std::log(x);
}
)";
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  sourceCodeHeader << std::endl << "/* This file was created by opendihu at " << StringUtility::timeToString(&tm)  //std::put_time(&tm, "%d/%m/%Y %H:%M:%S")
    << ".\n * It is designed for the FastMonodomainSolver and contains code for offloading to GPU.\n "
    << " */\n";

  VLOG(1) << "call defineHelperFunctions with helperFunctions: " << helperFunctions;

  // define helper functions
  sourceCodeHeader << defineHelperFunctions(helperFunctions, approximateExponentialFunction, false);
    
  // define main code that computes the rhs
  
  long int nInstancesToCompute = nInstancesToComputePerFiber * nFibersToCompute;
  std::string indent(10, ' ');
  
  sourceCodeMain << "\n" << indent << "// CellML define constants\n";
    
/*    << R"(  std::cout << "currentTime=" << currentTime << ", dt0D=" << dt0D << ", stimulate=" << stimulate << std::endl;)" << "\n" */
/*    << R"(  std::cout << "states[0]=" << states[0][0] << "," << states[0][1] << "," << states[0][2] << "," << states[0][3] << "," << std::endl;)" << "\n"
    << R"(  std::cout << "states[1]=" << states[1][0] << "," << states[1][1] << "," << states[1][2] << "," << states[1][3] << "," << std::endl;)" << "\n"
    << R"(  std::cout << "states[2]=" << states[2][0] << "," << states[2][1] << "," << states[2][2] << "," << states[2][3] << "," << std::endl;)" << "\n"
    << R"(  std::cout << "states[3]=" << states[3][0] << "," << states[3][1] << "," << states[3][2] << "," << states[3][3] << "," << std::endl;)" << "\n"*/
/*    << R"(  std::cout << "parameters[0]=" << parameters[0][0] << "," << parameters[0][1] << "," << parameters[0][2] << "," << parameters[0][3] << "," << std::endl;)" << "\n"
    << R"(  std::cout << "parameters[1]=" << parameters[1][0] << "," << parameters[1][1] << "," << parameters[1][2] << "," << parameters[1][3] << "," << std::endl;)" << "\n";*/

  // add assignments of constant values
  for (std::string constantAssignmentsLine : constantAssignments_)
  {
    constantAssignmentsLine = StringUtility::replaceAll(constantAssignmentsLine, "CONSTANTS[", "constant");
    constantAssignmentsLine = StringUtility::replaceAll(constantAssignmentsLine, "]", "");

    sourceCodeMain << indent << "const double " << constantAssignmentsLine << std::endl;
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
      bool first = true;
      
      codeExpression.visitLeafs([&sourceCodeLine,&isCommentedOut,nInstancesToCompute,this,&first](
        CellmlSourceCodeGeneratorVc::code_expression_t &expression, bool isFirstVariable)
      {
        switch(expression.type)
        {
          case code_expression_t::variableName:

            if (expression.code == "CONSTANTS")
            {
              // constants only exist once for all instances              
              if (first)
                sourceCodeLine << "const double ";
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
                //sourceCodeLine << "rate" << expression.arrayIndex;
                sourceCodeLine << "rates[" << expression.arrayIndex*nInstancesToCompute << "+instanceToComputeNo]";
              }
              else if (expression.code == "algebraics")
              {
                if (first)
                  sourceCodeLine << "const double ";
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
        first = false;
      });
      
      sourceCodeMain << indent << sourceCodeLine.str() << std::endl;
    }
  }
  sourceCodeMain << "\n"
    << "          // algebraic step\n"
    << "          // compute y* = y_n + dt*rhs(y_n), y_n = state, rhs(y_n) = rate, y* = intermediateState\n";

  //sourceCodeMain << indent << "double intermediateState0 = vmValues[instanceToComputeNo] + dt0D*rate0;\n";
  sourceCodeMain << indent << "states[0+instanceToComputeNo] = vmValues[instanceToComputeNo] + dt0D*rates[0+instanceToComputeNo];\n";

  for (int stateNo = 1; stateNo < this->nStates_; stateNo++)
  {
    sourceCodeMain << indent
       //<< "const double intermediateState" << stateNo << " = states[" << stateNo*nInstancesToCompute << "+instanceToComputeNo] + dt0D*rate" << stateNo << ";\n";
       << "states[" << stateNo*nInstancesToCompute << "+instanceToComputeNo] += dt0D*rates[" << stateNo*nInstancesToCompute << "+instanceToComputeNo];\n";
  }
  
  sourceCodeMain << "\n"
      << R"(
          // if stimulation, set value of Vm (state0)
          if (stimulateCurrentPoint)
          {
            //intermediateState0 = valueForStimulatedPoint;
             states[instanceToComputeNo] = valueForStimulatedPoint;
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
      bool first = true;
      
      codeExpression.visitLeafs([&sourceCodeLine,&isCommentedOut,nInstancesToCompute,&first,this](CellmlSourceCodeGeneratorVc::code_expression_t &expression, bool isFirstVariable)
      {
        switch(expression.type)
        {
          case code_expression_t::variableName:

            if (expression.code == "CONSTANTS")
            {
              // constants only exist once for all instances
              if (first)
                sourceCodeLine << "const double ";
              sourceCodeLine << "constant" << expression.arrayIndex;
            }
            else
            {
              // all other variables (states, rates, algebraics, parameters) exist for every instance
              if (expression.code == "states")
              {
                //sourceCodeLine << "intermediateState" << expression.arrayIndex;
                sourceCodeLine << "states[" << expression.arrayIndex*nInstancesToCompute << "+instanceToComputeNo]";
              }
              else if (expression.code == "rates")
              {
                //sourceCodeLine << "intermediateRate" << expression.arrayIndex;
                sourceCodeLine << "intermediateRates[" << expression.arrayIndex*nInstancesToCompute << "+instanceToComputeNo]";
              }
              else if (expression.code == "algebraics")
              {
                if (first)
                  sourceCodeLine << "const double ";
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
        first = false;
      });
      sourceCodeMain << indent << sourceCodeLine.str() << std::endl;
    }
  }

  sourceCodeMain << R"(
          // final step
          // y_n+1 = y_n + 0.5*[rhs(y_n) + rhs(y*)]
)";

  sourceCodeMain << indent << "vmValues[instanceToComputeNo] += 0.5*dt0D*("
    << "rates[0+instanceToComputeNo] " 
      << "+ intermediateRates[0+instanceToComputeNo]);\n";
  
  for (int stateNo = 1; stateNo < this->nStates_; stateNo++)
  {
    sourceCodeMain << indent << "states[" << stateNo*nInstancesToCompute << "+instanceToComputeNo] += 0.5*dt0D*(" 
      << "rates[" << stateNo*nInstancesToCompute << "+instanceToComputeNo] " 
      << "+ intermediateRates[" << stateNo*nInstancesToCompute << "+instanceToComputeNo]);\n";
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
  
  additionalCompileFlags_ = "-fopenmp -foffload=\"-O0 -g -v -lm\"";
  compilerCommand_ = CXX_COMPILER_COMMAND;
  sourceFileSuffix_ = ".cpp";
}
