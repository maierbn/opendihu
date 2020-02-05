#include "cellml/source_code_generator/source_code_generator.h"

#include <Python.h>  // has to be the first included header

#include "utility/string_utility.h"
#include "output_writer/generic.h"

#include <vector>
#include <iostream>
#include <Vc/Vc>
#include "easylogging++.h"

void CellMLSourceCodeGenerator::
generateSourceFileVc(std::string outputFilename, bool approximateExponentialFunction)
{
  //LOG(WARNING) << "In CellmlAdapter optimizationType \"vc\" is not yet implemented, fallback to \"simd\".";
  //generateSourceFileSimdAutovectorization(outputFilename);

  // replace pow and ?: functions

  std::set<std::string> helperFunctions;   //< functions found in the CellML code that need to be provided, usually the pow2, pow3, etc. helper functions for pow(..., 2), pow(...,3) etc.

  // loop over lines of CellML code
  for (code_expression_t &codeExpression : cellMLCode_.lines)
  {
    VLOG(1) << "line: " << codeExpression.getString();

    if (codeExpression.type != code_expression_t::commented_out)
    {
      // loop over all nodes in the syntax tree of this line
      codeExpression.visitNodes([&helperFunctions](CellMLSourceCodeGenerator::code_expression_t &expression)
      {
        // we look for occurences of functions and ternary operators, these can be detected in a tree node
        if (expression.type == code_expression_t::tree)
        {

          std::stringstream a;
          for (int i = 0; i < expression.treeChildren.size(); i++)
          {
            a << "[" << expression.treeChildren[i].code << "] ";
          }
          VLOG(1) << "check expression: tree with " << expression.treeChildren.size() << " children: " << a.str();

          // loop over children of the tree node
          for (int i = 0; i < expression.treeChildren.size(); i++)
          {
            code_expression_t &innerExpression = expression.treeChildren[i];

            if (innerExpression.type == code_expression_t::otherCode)
            {
              // if innerExpression contains start of "pow" function, find closing ")", extract exponent (e.g. pow(..., 3.000) -> 3) and replace by pow3 function
              if (innerExpression.code.find("pow") != std::string::npos
                  && innerExpression.code.find("pow") == innerExpression.code.length()-3)
              {
                int exponent = 0;

                VLOG(1) << "found pow at i=" << i << " of [" << innerExpression.code << "]";

                // now the next expressions should be "(" "<anything>,exponent" ")"
                assert(expression.treeChildren.size() > i + 1);
                assert(expression.treeChildren[i+1].type == code_expression_t::tree);
                assert(expression.treeChildren[i+1].treeChildren.size() == 3);

                assert(expression.treeChildren[i+1].treeChildren[0].type == code_expression_t::otherCode);
                assert(expression.treeChildren[i+1].treeChildren[0].code == "(");
                assert(expression.treeChildren[i+1].treeChildren[1].type == code_expression_t::tree);
                assert(expression.treeChildren[i+1].treeChildren[2].type == code_expression_t::otherCode);
                assert(expression.treeChildren[i+1].treeChildren[2].code == ")");

                code_expression_t &expressionExponent = expression.treeChildren[i+1].treeChildren[1].treeChildren.back();

                if (expressionExponent.type == code_expression_t::otherCode)
                {
                  std::string code = expressionExponent.code;
                  std::size_t posComma = code.find(",");

                  assert(posComma != std::string::npos);

                  std::string codeExponent = code.substr(posComma+1);
                  StringUtility::trim(codeExponent);
                  exponent = atoi(codeExponent.c_str());

                  // remove ", exponent" from code
                  expressionExponent.code = code.substr(0, posComma);
                }
                else if (expressionExponent.type == code_expression_t::tree)
                {
                  std::string code = expressionExponent.treeChildren.back().code;
                  std::size_t posComma = code.find(",");

                  assert(posComma != std::string::npos);

                  std::string codeExponent = code.substr(posComma+1);
                  StringUtility::trim(codeExponent);
                  exponent = atoi(codeExponent.c_str());

                  // remove ", exponent" from code
                  expressionExponent.treeChildren.back().code = code.substr(0, posComma);
                }

                // compose name of helper function: pow<exponent>
                std::stringstream s;
                if (exponent > 0)
                {
                  s << "pow" << exponent;
                }
                else
                {
                  // a negative exponent, e.g. "pow(x,-2)" yields "powReciprocal2(x)"
                  s << "powReciprocal" << -exponent;
                }
                std::string helperFunction = s.str();
                helperFunctions.insert(helperFunction);

                // replace "pow" by e.g. "pow3" for exponent 3
                std::size_t posPow = innerExpression.code.find("pow");
                innerExpression.code = innerExpression.code.substr(0, posPow) + helperFunction;
              }
              else if (innerExpression.code == "?")
              {
                assert(expression.treeChildren.size() == 5);  //<condition> "?" <branch0> ":" <branch1>

                code_expression_t iifFunction;
                iifFunction.type = code_expression_t::tree;
                iifFunction.treeChildren.resize(8);

                iifFunction.treeChildren[0].type = code_expression_t::otherCode;
                iifFunction.treeChildren[0].code = "Vc::iif";

                iifFunction.treeChildren[1].type = code_expression_t::otherCode;
                iifFunction.treeChildren[1].code = "(";

                iifFunction.treeChildren[2] = expression.treeChildren[0];

                iifFunction.treeChildren[3].type = code_expression_t::otherCode;
                iifFunction.treeChildren[3].code = ",";

                iifFunction.treeChildren[4] = expression.treeChildren[2];

                iifFunction.treeChildren[5].type = code_expression_t::otherCode;
                iifFunction.treeChildren[5].code = ",";

                iifFunction.treeChildren[6] = expression.treeChildren[4];

                iifFunction.treeChildren[7].type = code_expression_t::otherCode;
                iifFunction.treeChildren[7].code = ")";

                // replace 0.00000 and 1.00000 by double_v(Vc::One) and double_v(Vc::Zero)

                // only do replace if
                bool enableReplace = true;

                // for k in {4,6}
                for (int k = 4; k <= 6; k +=2)
                {
                  bool expressionHasVariables = false;
                  iifFunction.treeChildren[k].visitLeafs([&expressionHasVariables](code_expression_t &expression, bool isFirstVariable)
                  {
                    if (expression.type == code_expression_t::variableName)
                    {
                      if (expression.code != "CONSTANTS")
                        expressionHasVariables = true;
                    }
                  });

                  if (iifFunction.treeChildren[k].type != code_expression_t::otherCode && !expressionHasVariables)
                    enableReplace = false;
                }

                if (enableReplace)
                {
                  // for k in {4,6}
                  for (int k = 4; k <= 6; k +=2)
                  {
                    if (iifFunction.treeChildren[k].type == code_expression_t::otherCode)
                    {
                      std::string code = iifFunction.treeChildren[k].code;
                      StringUtility::trim(code);

                      if (code == "0.00000")
                      {
                        iifFunction.treeChildren[k].code = " Vc::double_v(Vc::Zero)";
                      }
                      else if (code == "1.00000")
                      {
                        iifFunction.treeChildren[k].code = " Vc::double_v(Vc::One)";
                      }
                      else if (atof(code.c_str()) != 0.0)
                      {
                        iifFunction.treeChildren[k].code = " Vc::double_v(Vc::One)*" + iifFunction.treeChildren[k].code;
                      }
                    }
                  }
                }

                expression = iifFunction;
                break;
              }
            }
          }
        }
        else if (expression.type == code_expression_t::otherCode)
        {
          if (expression.code.find("exp") != std::string::npos
              && expression.code.find("exp") == expression.code.length()-3)
          {
            // replace "exp(" by "exponential("
            expression.code = StringUtility::replace(expression.code, "exp", "exponential");
            helperFunctions.insert("exponential");
          }
        }
      });
    }
  }

  VLOG(1) << "helperFunctions: " << helperFunctions;

  std::stringstream simdSource;
  simdSource << "#include <math.h>" << std::endl
    << "#include <Vc/Vc>" << std::endl
    << cellMLCode_.header << std::endl
    << "using Vc::double_v; " << std::endl;

  // define helper functions

  // add helper functions for helper functions (e.g. pow4 needs pow2)
  int previousSize = -1;
  while (helperFunctions.size() != previousSize)
  {
    previousSize = helperFunctions.size();

    for (std::set<std::string>::iterator iter = helperFunctions.begin(); iter != helperFunctions.end(); iter++)
    {
      std::string functionName = *iter;
      if (functionName.find("pow") != std::string::npos)
      {
        int exponent = atoi(functionName.substr(3).c_str());

        if (functionName.find("powReciprocal") != std::string::npos)
          exponent = atoi(functionName.substr(std::string("powReciprocal").length()).c_str());

        int exponent0 = int(exponent/2);
        int otherExponent = exponent - exponent0;

        std::stringstream requiredFunction;
        requiredFunction << "pow" << exponent0;

        if (exponent0 != 0 && exponent0 != 1)
          helperFunctions.insert(requiredFunction.str());

        requiredFunction.str("");
        requiredFunction << "pow" << otherExponent;
        if (otherExponent != 1)
          helperFunctions.insert(requiredFunction.str());
      }
    }
  }

  // generate declarations
  simdSource << "\n// helper functions\n";
  for (std::set<std::string>::iterator iter = helperFunctions.begin(); iter != helperFunctions.end(); iter++)
  {
    std::string functionName = *iter;
    simdSource << "Vc::double_v " << functionName << "(Vc::double_v x);" << std::endl;
  }

  // define exp function if needed
  if (helperFunctions.find("exponential") != helperFunctions.end())
  {
    if (approximateExponentialFunction)
    {
      simdSource << R"(
Vc::double_v exponential(Vc::double_v x)
{
  //return Vc::exp(x);
  // it was determined the x is always in the range [-12,+12]

  // exp(x) = lim n→∞ (1 + x/n)^n, we set n=1024
  x = 1.0 + x / 1024.;
  for (int i = 0; i < 10; i++)
  {
    x *= x;
  }
  return x;

  // relative error of this implementation:
  // x    rel error
  // 0    0
  // 1    0.00048784455634225593
  // 3    0.0043763626896140342
  // 5    0.012093715791500804
  // 9    0.038557535762274039
  // 12   0.067389808619653505
}
)";
    }
    else
    {
      simdSource << R"(
Vc::double_v exponential(Vc::double_v x)
{
  return Vc::exp(x);
}
)";
    }
  }

  // generate other helper functions
  for (std::set<std::string>::iterator iter = helperFunctions.begin(); iter != helperFunctions.end(); iter++)
  {
    std::string functionName = *iter;

    // exp function is already handled
    if (functionName == "exponential")
      continue;

    if (functionName.find("pow") != std::string::npos)
    {
      int exponent = atoi(functionName.substr(3).c_str());
      if (functionName.find("powReciprocal") != std::string::npos)
      {
        exponent = -atoi(functionName.substr(std::string("powReciprocal").length()).c_str());
      }
      if (exponent == 2)
      {
      simdSource << R"(
Vc::double_v pow2(Vc::double_v x)
{
  return x*x;
}
)";
      }
      else
      {
        int exponent0 = int(fabs(exponent)/2);
        int otherExponent = fabs(exponent) - exponent0;
        simdSource << "Vc::double_v pow" << (exponent < 0? "Reciprocal" : "")
          << fabs(exponent) << "(Vc::double_v x)" << std::endl
          << "{" << std::endl
          << "  return ";
        if (exponent < 0)
          simdSource << "1./(";

        // if exponent == 1 => exponent0 == 0
        if (exponent0 == 0)
        {
          simdSource << "x";
        }
        else
        {
          if (exponent0 == 1)
          {
            simdSource << "x*(";
          }
          else
          {
            simdSource << "pow" << exponent0 << "(";
          }
          simdSource << "pow" << otherExponent << "(x))";
        }

        if (exponent < 0)
          simdSource << ")";
        simdSource << ";" << std::endl;

        simdSource << "}" << std::endl << std::endl;
      }
    }
  }


  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  simdSource << std::endl << "/* This function was created by opendihu at " << StringUtility::timeToString(&tm)  //std::put_time(&tm, "%d/%m/%Y %H:%M:%S")
    << ".\n * It is designed for " << this->nInstances_ << " instances of the CellML problem.\n "
    << "The \"optimizationType\" is \"vc\". (Other options are \"simd\" and \"openmp\".) */" << std::endl
    << "extern \"C\"" << std::endl
    << "void computeCellMLRightHandSide("
    << "void *context, double t, double *states, double *rates, double *intermediates, double *parameters)" << std::endl << "{" << std::endl;

  simdSource << "  double VOI = t;   /* current simulation time */" << std::endl;
  simdSource << std::endl << "  /* define constants */" << std::endl
    << "  double CONSTANTS[" << this->nConstants_ << "];" << std::endl;

  // add assignments of constant values
  for (std::string constantAssignmentsLine : constantAssignments_)
  {
    simdSource << "  " << constantAssignmentsLine << std::endl;
  }

  // add declaration of algebraic variables
  simdSource << std::endl;
  const int nVcVectors = (int)(ceil((double)this->nInstances_ / Vc::double_v::Size));
  const int nParametersPerInstance = this->parameters_.size() / this->nInstances_;
  const int vcVectorSize = (int)Vc::double_v::Size;

  simdSource << std::endl
    << "  const int nInstances = " << this->nInstances_ << ";\n"
    << "  const int nStates = " << this->nStates_ << ";\n"
    << "  const int nIntermediates = " << this->nIntermediates_ << ";\n"
    << "  const int nParametersPerInstance = " << nParametersPerInstance << ";\n"
    << "\n"
    << "  const int vcVectorSize = " << (int)Vc::double_v::Size << ";\n"
    << "  const int nVcVectors = " << nVcVectors << ";  // ceil(" << this->nInstances_ << " instances / VcSize " << Vc::double_v::Size << ")" << std::endl
    << "\n"
    << "  if (vcVectorSize != " << vcVectorSize << ")\n"
    << "  {\n"
    << "    std::cout << \"Error! SIMD register length mismatches between opendihu (" << vcVectorSize
      << ") and compiled code \" << vcVectorSize << \"!\" << std::endl;\n"
    << "    return;\n"
    << "  }\n"
    << "  Vc::double_v statesVc[" << nStates_*nVcVectors << "];  // " << this->nStates_ << " states * " << nVcVectors << " vectors" << std::endl
    << "  Vc::double_v ratesVc[" << nStates_*nVcVectors << "];   // " << this->nStates_ << " rates  * " << nVcVectors << " vectors" << std::endl
    << "  Vc::double_v intermediatesVc[" << nIntermediates_*nVcVectors << "];  // " << this->nIntermediates_ << " intermediates  * " << nVcVectors << " vectors" << std::endl
    << "  Vc::double_v parametersVc[" << nParametersPerInstance*nVcVectors << "];  // " << nParametersPerInstance << " parameters  * " << nVcVectors << " vectors" << std::endl
    << "\n"
    << "  // fill input vectors of states and parameters\n"
    << "  for (int stateNo = 0; stateNo < nStates; stateNo++)\n"
    << "    for (int i = 0; i < nVcVectors; i++)  // Vc vector no\n"
    << "      for (int k = 0; k < vcVectorSize; k++)  // entry no in Vc vector \n"
    << "      {\n"
    << "        if (i*" << vcVectorSize << "+k >= " << nInstances_ << ")\n"
    << "          continue;\n"
    << "        int index = stateNo*" << nInstances_ << " + i*" << vcVectorSize << "+k;\n"
    << "        statesVc[stateNo*" << nVcVectors << " + i][k] = states[index];\n"
    << "      }\n"
    << "\n"
    << "  for (int parameterNo = 0; parameterNo < nParametersPerInstance; parameterNo++)\n"
    << "    for (int i = 0; i < nVcVectors; i++)  // Vc vector no\n"
    << "      for (int k = 0; k < vcVectorSize; k++)  // entry no in Vc vector \n"
    << "      {\n"
    << "        if (i*" << vcVectorSize << "+k >= " << nInstances_ << ")\n"
    << "          continue;\n"
    << "        int index = parameterNo*" << nInstances_ << " + i*" << vcVectorSize << "+k;\n"
    << "        parametersVc[parameterNo*" << nVcVectors << " + i][k] = parameters[index];\n"
    << "      }\n"
    << std::endl
    << "  for (int i = 0; i < nVcVectors; i++)" << std::endl
    << "  {" << std::endl;

  std::stringstream debug;

  // loop over lines of cellml code
  for (code_expression_t &codeExpression : cellMLCode_.lines)
  {
    if (codeExpression.type != code_expression_t::commented_out)
    {
      simdSource << "    ";
      debug << "    ";

      codeExpression.visitLeafs([&simdSource,&nVcVectors,&debug,this](CellMLSourceCodeGenerator::code_expression_t &expression, bool isFirstVariable)
      {
        switch(expression.type)
        {
          case code_expression_t::variableName:

            if (expression.code == "CONSTANTS")
            {
              // constants only exist once for all instances
              simdSource << expression.code << "[" << expression.arrayIndex<< "]";
              debug << expression.code << "[" << expression.arrayIndex<< "]";
            }
            else
            {
              // all other variables (states, rates, intermediates, parameters) exist for every instance
              if (expression.code == "states")
              {
                simdSource << "statesVc[" << expression.arrayIndex * nVcVectors << "+i]";
                debug << "states[" << expression.arrayIndex * this->nInstances_ << "+i]";
              }
              else if (expression.code == "rates")
              {
                simdSource << "ratesVc[" << expression.arrayIndex * nVcVectors << "+i]";
                debug << "rates[" << expression.arrayIndex * this->nInstances_ << "+i]";
              }
              else if (expression.code == "intermediates")
              {
                simdSource << "intermediatesVc[" << expression.arrayIndex * nVcVectors << "+i]";
                debug << "intermediates[" << expression.arrayIndex * this->nInstances_ << "+i]";
              }
              else if (expression.code == "parameters")
              {
                simdSource << "parametersVc[" << expression.arrayIndex * nVcVectors << "+i]";
                debug << "parameters[" << expression.arrayIndex * this->nInstances_ << "+i]";
              }
              else
              {
                LOG(FATAL) << "unhandled variable type \"" << expression.code << "\".";
                simdSource << expression.code << "[" << expression.arrayIndex << " * nInstances + i]";
              }
            }
            break;

          case code_expression_t::otherCode:
            simdSource << expression.code;
            debug << expression.code;
            break;

          case code_expression_t::commented_out:
            simdSource << "  // (not assigning to a parameter) " << expression.code;
            debug << "  // (not assigning to a parameter) " << expression.code;
            break;

          default:
            break;
        }
      });
      simdSource << std::endl;
      debug << std::endl;
    }
  }
  simdSource << std::endl;
  //  << "/*\n" << debug.str() << "*/\n\n";

  simdSource << "  }" << std::endl << std::endl
    << "  // store computed values back to pointers\n"
    << "  for (int rateNo = 0; rateNo < nStates; rateNo++)\n"
    << "    for (int i = 0; i < nVcVectors; i++)  // Vc vector no\n"
    << "      for (int k = 0; k < " << vcVectorSize << "; k++)  // entry no in Vc vector \n"
    << "      {\n"
    << "        if (i*" << vcVectorSize << "+k >= " << nInstances_ << ")\n"
    << "          continue;\n"
    << "        int index = rateNo*" << nInstances_ << " + i*" << vcVectorSize << "+k;\n"
    << "        rates[index] = ratesVc[rateNo*" << nVcVectors << " + i][k];\n"
    << "      }\n"
    << "\n"
    << "  for (int intermediateNo = 0; intermediateNo < nIntermediates; intermediateNo++)\n"
    << "    for (int i = 0; i < nVcVectors; i++)  // Vc vector no\n"
    << "      for (int k = 0; k < " << vcVectorSize << "; k++)  // entry no in Vc vector \n"
    << "      {\n"
    << "        if (i*" << vcVectorSize << "+k >= " << nInstances_ << ")\n"
    << "          continue;\n"
    << "        int index = intermediateNo*" << nInstances_ << "+ i*" << vcVectorSize << "+k;\n"
    << "        intermediates[index] = intermediatesVc[intermediateNo*" << nVcVectors << " + i][k];\n"
    << "      }\n"
    << "\n";


  // add footer
  simdSource << cellMLCode_.footer << std::endl;

  // write out source file
  std::ofstream simdSourceFile;
  OutputWriter::Generic::openFile(simdSourceFile, outputFilename);
  if (!simdSourceFile.is_open())
  {
    LOG(FATAL) << "Could not write to file \"" << outputFilename << "\".";
  }
  else
  {
    std::string fileContents = simdSource.str();
    simdSourceFile << fileContents;
    simdSourceFile.close();
  }

  std::stringstream s;
  s << "-lVc -I\"" << OPENDIHU_HOME << "/dependencies/vc/install/include\" "
    << "-L\"" << OPENDIHU_HOME << "/dependencies/vc/install/lib\" ";
  additionalCompileFlags_ = s.str();
  compilerCommand_ = CXX_COMPILER_COMMAND;
  sourceFileSuffix_ = ".cpp";
}
