#include "cellml/source_code_generator/02_generator_openmp.h"

#include <Python.h>  // has to be the first included header

#include "utility/string_utility.h"
#include "output_writer/generic.h"

#include <vector>
#include <iostream>
#include <Vc/Vc>
#include "easylogging++.h"

void CellmlSourceCodeGeneratorOpenMp::
generateSourceFileOpenMP(std::string outputFilename, int maximumNumberOfThreads)
{
  std::stringstream sourceCode;
  sourceCode << "#include <math.h>" << std::endl
    << "#include <omp.h>" << std::endl
    << cellMLCode_.header << std::endl;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  sourceCode << std::endl << "/* This function was created by opendihu at " << StringUtility::timeToString(&tm)  //std::put_time(&tm, "%d/%m/%Y %H:%M:%S")
    << ".\n * It is designed for " << this->nInstances_ << " instances of the CellML problem.\n "
    << " * The \"optimizationType\" is \"openmp\". (Other options are \"vc\" and \"simd\".) */" << std::endl
    << "void computeCellMLRightHandSide("
    << "void *context, double t, double *states, double *rates, double *intermediates, double *parameters)" << std::endl << "{" << std::endl;

  if (maximumNumberOfThreads > 0)
  {
    sourceCode << "  omp_set_num_threads(" << maximumNumberOfThreads << ");\n";
  }

  sourceCode << "  double VOI = t;   /* current simulation time */" << std::endl;
  sourceCode << std::endl << "  /* define constants */" << std::endl
    << "  double CONSTANTS[" << this->nConstants_ << "];" << std::endl;

  // add assignments of constant values
  for (std::string constantAssignmentsLine : constantAssignments_)
  {
    sourceCode << "  " << constantAssignmentsLine << std::endl;
  }

  sourceCode << std::endl
    << "  #pragma omp parallel for" << std::endl
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
            // all other variables (states, rates, intermediates, parameters) exist for every instance
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

  additionalCompileFlags_ = "-fopenmp";
  compilerCommand_ = C_COMPILER_COMMAND;
  sourceFileSuffix_ = ".c";
}
