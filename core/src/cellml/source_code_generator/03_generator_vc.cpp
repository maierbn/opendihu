#include "cellml/source_code_generator/03_generator_vc.h"

#include <Python.h>  // has to be the first included header

#include "utility/string_utility.h"
#include "output_writer/generic.h"

#include <vector>
#include <iostream>
#include <Vc/Vc>
#include "easylogging++.h"

void CellmlSourceCodeGeneratorVc::preprocessCode(std::set<std::string> &helperFunctions)
{
  if (preprocessingDone_)
    return;

  // loop over lines of CellML code
  for (code_expression_t &codeExpression : cellMLCode_.lines)
  {
    VLOG(1) << "line: " << codeExpression.getString();

    if (codeExpression.type != code_expression_t::commented_out)
    {
      // loop over all nodes in the syntax tree of this line
      codeExpression.visitNodes([&helperFunctions](CellmlSourceCodeGeneratorVc::code_expression_t &expression)
      {
        // we look for occurences of functions and ternary operators, these can be detected in a tree node
        if (expression.type == code_expression_t::tree)
        {
          if (VLOG_IS_ON(1))
          {
            std::stringstream a;
            for (int i = 0; i < expression.treeChildren.size(); i++)
            {
              a << "[" << expression.treeChildren[i].code << "] ";
            }
            VLOG(1) << "check expression: tree with " << expression.treeChildren.size() << " children: " << a.str();
          }

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
                bool isIntegerExponent = false;

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
                  // exponent can also be "- 1.0000", so remove all whitespace in the inner
                  codeExponent.erase(std::remove(codeExponent.begin(), codeExponent.end(), ' '), codeExponent.end());
                  
                  isIntegerExponent = codeExponent.find_first_not_of("0123456789+-") == std::string::npos;
                  exponent = atoi(codeExponent.c_str());
                  if (exponent == 0)
                    isIntegerExponent = false;
                    
                  if (isIntegerExponent)
                  {
                    // remove ", exponent" from code
                    expressionExponent.code = code.substr(0, posComma);
                  }
                }
                else if (expressionExponent.type == code_expression_t::tree)
                {
                  std::string code = expressionExponent.treeChildren.back().code;
                  std::size_t posComma = code.find(",");

                  assert(posComma != std::string::npos);

                  std::string codeExponent = code.substr(posComma+1);
                  StringUtility::trim(codeExponent);
                  // exponent can also be "- 1.0000", so remove all whitespace in the inner
                  codeExponent.erase(std::remove(codeExponent.begin(), codeExponent.end(), ' '), codeExponent.end());
                  
                  isIntegerExponent = codeExponent.find_first_not_of("0123456789+-") == std::string::npos;
                  exponent = atoi(codeExponent.c_str());
                  if (exponent == 0)
                    isIntegerExponent = false;

                  if (isIntegerExponent)
                  {
                    // remove ", exponent" from code
                    expressionExponent.treeChildren.back().code = code.substr(0, posComma);
                  }
                }
                
                // if the exponent of the pow function is just a single integer number
                if (isIntegerExponent)
                {
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
                else 
                {
                  helperFunctions.insert("pow");
                }
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

                // only do replace if not both `then` and `else` contain variables
                bool enableReplace = false;

                // for k in {4,6}
                for (int k = 4; k <= 6; k +=2)
                {
                  // check if the code contains variables, then it should not be replaced because it is already vectorized
                  bool expressionHasVariables = false;
                  iifFunction.treeChildren[k].visitLeafs([&expressionHasVariables](code_expression_t &expression, bool isFirstVariable)
                  {
                    if (expression.type == code_expression_t::variableName)
                    {
                      if (expression.code != "CONSTANTS")   // a constant does not count as a variable here
                        expressionHasVariables = true;
                    }
                  });

                  if (iifFunction.treeChildren[k].type == code_expression_t::otherCode 
                      || !expressionHasVariables)
                    enableReplace = true;      // replacing is necessary because in the current slot there are no variables
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
                      // code can also be "- 1.0000", so remove all whitespace in the inner
                      code.erase(std::remove(code.begin(), code.end(), ' '), code.end());

                      if (code == "0.00000")
                      {
                        iifFunction.treeChildren[k].code = " Vc::double_v(Vc::Zero)";
                      }
                      else if (code == "1.00000")
                      {
                        iifFunction.treeChildren[k].code = " Vc::double_v(Vc::One)";
                      }
                      else if (atof(code.c_str()) != 0.0)   // if there is numeric code, not just '('s
                      {
                        iifFunction.treeChildren[k].code = " (Vc::double_v(Vc::One)*(" + iifFunction.treeChildren[k].code + "))";
                      }
                    }
                    else if (iifFunction.treeChildren[k].type == code_expression_t::tree)
                    {
                      iifFunction.treeChildren[k].visitLeafs([](code_expression_t &expression, bool isFirstVariable)
                      {
                        if (expression.type == code_expression_t::variableName)
                        {
                          if (expression.code == "CONSTANTS")
                          {
                            // replace a single constant
                            expression.type = code_expression_t::tree;
                            expression.treeChildren.resize(3);

                            expression.treeChildren[0].type = code_expression_t::otherCode;
                            expression.treeChildren[0].code = "(Vc::double_v(Vc::One)*";

                            expression.treeChildren[1].type = code_expression_t::variableName;
                            expression.treeChildren[1].code = "CONSTANTS";
                            expression.treeChildren[1].arrayIndex = expression.arrayIndex;
                            
                            expression.treeChildren[2].type = code_expression_t::otherCode;
                            expression.treeChildren[2].code = ")";
                          }
                        }
                      });
                    }
                    else if (iifFunction.treeChildren[k].type == code_expression_t::variableName
                             && iifFunction.treeChildren[k].code == "CONSTANTS")
                    {
                      // replace a single constant
                      iifFunction.treeChildren[k].type = code_expression_t::tree;
                      iifFunction.treeChildren[k].treeChildren.resize(3);

                      iifFunction.treeChildren[k].treeChildren[0].type = code_expression_t::otherCode;
                      iifFunction.treeChildren[k].treeChildren[0].code = "(Vc::double_v(Vc::One)*";

                      iifFunction.treeChildren[k].treeChildren[1].type = code_expression_t::variableName;
                      iifFunction.treeChildren[k].treeChildren[1].code = "CONSTANTS";
                      iifFunction.treeChildren[k].treeChildren[1].arrayIndex = iifFunction.treeChildren[k].arrayIndex;
                      
                      iifFunction.treeChildren[k].treeChildren[2].type = code_expression_t::otherCode;
                      iifFunction.treeChildren[k].treeChildren[2].code = ")";
                    }
                  }
                }

                expression = iifFunction;
                break;
              }
              else if (innerExpression.code == "fabs")
              {
                innerExpression.code = "Vc::abs";
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

  preprocessingDone_ = true;
}

std::string CellmlSourceCodeGeneratorVc::defineHelperFunctions(std::set<std::string> &helperFunctions, bool approximateExponentialFunction)
{
  if (!helperFunctionsCode_.empty())
    return helperFunctionsCode_;

  std::stringstream sourceCode;

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
        if (exponent == 0)
        {
          helperFunctions.insert("exponential");
          continue;
        }

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

  VLOG(1) << "after adding all necessary helperFunctions: " << helperFunctions;
  
  // generate declarations
  sourceCode << "\n// helper functions\n";
  for (std::set<std::string>::iterator iter = helperFunctions.begin(); iter != helperFunctions.end(); iter++)
  {
    std::string functionName = *iter;
    if (functionName == "pow")
    {
      sourceCode << "Vc::double_v pow(Vc::double_v basis, Vc::double_v exponent);" << std::endl;
      sourceCode << "Vc::double_v pow(Vc::double_v basis, double exponent);" << std::endl;
    }
    else
    {
      sourceCode << "Vc::double_v " << functionName << "(Vc::double_v x);" << std::endl;
    }
  }

  // define exp function if needed
  if (helperFunctions.find("exponential") != helperFunctions.end())
  {
    if (approximateExponentialFunction)
    {
      sourceCode << R"(
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
      sourceCode << R"(
Vc::double_v exponential(Vc::double_v x)
{
  return Vc::exp(x);
}
)";
    }
  }
  
  // define pow function if needed
  if (helperFunctions.find("pow") != helperFunctions.end())
  {
    sourceCode << R"(
Vc::double_v pow(Vc::double_v basis, Vc::double_v exponent)
{
  // Note, there is no pow function defined by Vc.
  return exponential(Vc::log(basis)*exponent);
}

Vc::double_v pow(Vc::double_v basis, double exponent)
{
  return exponential(Vc::log(basis)*exponent);
}

)";
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
      if (exponent != 0)
      {
        if (functionName.find("powReciprocal") != std::string::npos)
        {
          exponent = -atoi(functionName.substr(std::string("powReciprocal").length()).c_str());
        }
        if (exponent == 2)
        {
        sourceCode << R"(
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
          sourceCode << "Vc::double_v pow" << (exponent < 0? "Reciprocal" : "")
            << fabs(exponent) << "(Vc::double_v x)" << std::endl
            << "{" << std::endl
            << "  return ";
          if (exponent < 0)
            sourceCode << "1./(";

          // if exponent == 1 => exponent0 == 0
          if (exponent0 == 0)
          {
            sourceCode << "x";
          }
          else
          {
            if (exponent0 == 1)
            {
              sourceCode << "x*(";
            }
            else
            {
              sourceCode << "pow" << exponent0 << "(";
            }
            sourceCode << "pow" << otherExponent << "(x))";
          }

          if (exponent < 0)
            sourceCode << ")";
          sourceCode << ";" << std::endl;

          sourceCode << "}" << std::endl << std::endl;
        }
      }
    }
  }

  helperFunctionsCode_ = sourceCode.str();
  VLOG(1) << "Helper functions produce the following code: " << helperFunctionsCode_;

  return helperFunctionsCode_;
}

void CellmlSourceCodeGeneratorVc::
generateSourceFileVc(std::string outputFilename, bool approximateExponentialFunction)
{
  std::set<std::string> helperFunctions;   //< functions found in the CellML code that need to be provided, usually the pow2, pow3, etc. helper functions for pow(..., 2), pow(...,3) etc.

  // replace pow and ?: functions
  preprocessCode(helperFunctions);

  VLOG(1) << "helperFunctions: " << helperFunctions;

  std::stringstream sourceCode;
  sourceCode << "#include <math.h>" << std::endl
    << "#include <Vc/Vc>" << std::endl
    << cellMLCode_.header << std::endl
    << "using Vc::double_v; " << std::endl;

  // define helper functions
  sourceCode << defineHelperFunctions(helperFunctions, approximateExponentialFunction);

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  sourceCode << std::endl << "// This function was created by opendihu at " << StringUtility::timeToString(&tm)  //std::put_time(&tm, "%d/%m/%Y %H:%M:%S")
    << ".\n// It is designed for " << this->nInstances_ << " instances of the CellML problem.\n"
    << "// The \"optimizationType\" is \"vc\". (Other options are \"simd\" and \"openmp\".)" << std::endl;
    
  if (!parametersUsedAsAlgebraic_.empty())
  {
    sourceCode << "// " << parametersUsedAsAlgebraic_.size() << " algebraic" << (parametersUsedAsAlgebraic_.size()==1? " is": "s are") << " replaced by parameters: ";
    bool first = true;
    for (int parameterUsedAsAlgebraic : parametersUsedAsAlgebraic_)
    {
      if (!first)
        sourceCode << ", ";
      sourceCode << algebraicNames_[parameterUsedAsAlgebraic];
      first = false;
    }
    sourceCode << std::endl;
  }
  if (!parametersUsedAsConstant_.empty())
  {
    sourceCode << "// " << parametersUsedAsConstant_.size() << " constant" << (parametersUsedAsConstant_.size()==1? " is": "s are") << " replaced by parameters: ";
    bool first = true;
    for (int parameterUsedAsConstant : parametersUsedAsConstant_)
    {
      if (!first)
        sourceCode << ", ";
      sourceCode << constantNames_[parameterUsedAsConstant];
      first = false;
    }
    sourceCode << std::endl;
  }
  sourceCode  << std::endl;
  std::vector<int> parametersUsedAsConstant_;  //< explicitely defined parameters that will be copied to constants, this vector contains the indices of the constants

    
  sourceCode << "#ifdef __cplusplus\n" << "extern \"C\"\n" << "#endif\n" << std::endl
    << "void computeCellMLRightHandSide("
    << "void *context, double t, double *states, double *rates, double *algebraics, double *parameters)" << std::endl 
    << "{" << std::endl
    << "  // assert that Vc::double_v::Size is the same as in opendihu, otherwise there will be problems\n"
    << "  if (Vc::double_v::Size != " << Vc::double_v::Size << ")\n"
    << "  {\n"
    << "    std::cout << \"Fatal error in compiled library of source file \\\"" << outputFilename << "\\\", size of SIMD register in "
    << "compiled code (\" << Vc::double_v::Size << \") does not match opendihu code (" << Vc::double_v::Size << ").\" << std::endl;\n"
    << "    std::cout << \"Delete library such that it will be regenerated with the correct compile options!\" << std::endl;\n"
    << "    exit(1);\n"
    << "  }\n\n";

  sourceCode << "  double VOI = t;   /* current simulation time */" << std::endl;
  sourceCode << std::endl << "  /* define constants */" << std::endl
    << "  double CONSTANTS[" << this->nConstants_ << "];" << std::endl;

  // add assignments of constant values
  for (std::string constantAssignmentsLine : constantAssignments_)
  {
    sourceCode << "  " << constantAssignmentsLine << std::endl;
  }

  // add declaration of algebraic variables
  sourceCode << std::endl;
  const int nVcVectors = (int)(ceil((double)this->nInstances_ / Vc::double_v::Size));
  const int nParametersPerInstance = this->nAlgebraics_;

  sourceCode << std::endl
    << "  const int nInstances = " << this->nInstances_ << ";\n"
    << "  const int nStates = " << this->nStates_ << ";\n"
    << "  const int nAlgebraics = " << this->nAlgebraics_ << ";\n"
    << "  const int nParametersPerInstance = " << nParametersPerInstance << ";\n"
    << "  const int nVcVectors = " << nVcVectors << ";  // ceil(" << this->nInstances_ << " instances / VcSize " << Vc::double_v::Size << ")" << std::endl
    << "  Vc::double_v statesVc[nStates*nVcVectors];  // " << this->nStates_ << " states * " << nVcVectors << " vectors" << std::endl
    << "  Vc::double_v ratesVc[nStates*nVcVectors];   // " << this->nStates_ << " rates  * " << nVcVectors << " vectors" << std::endl
    << "  Vc::double_v algebraicsVc[nAlgebraics*nVcVectors];  // " << this->nAlgebraics_ << " algebraics  * " << nVcVectors << " vectors" << std::endl
    << "  Vc::double_v parametersVc[nParametersPerInstance*nVcVectors];  // " << nParametersPerInstance << " parameters  * " << nVcVectors << " vectors" << std::endl
    << "\n"
    << "  // fill input vectors of states and parameters\n"
    << "  for (int stateNo = 0; stateNo < nStates; stateNo++)\n"
    << "    for (int i = 0; i < nVcVectors; i++)  // Vc vector no\n"
    << "      for (int k = 0; k < Vc::double_v::Size; k++)  // entry no in Vc vector \n"
    << "        statesVc[stateNo*nVcVectors + i][k] = states[std::min(stateNo*nInstances + i*(int)Vc::double_v::Size+k, nStates*nInstances-1)];\n"
    << "\n"
    << "  for (int parameterNo = 0; parameterNo < nParametersPerInstance; parameterNo++)\n"
    << "    for (int i = 0; i < nVcVectors; i++)  // Vc vector no\n"
    << "      for (int k = 0; k < Vc::double_v::Size; k++)  // entry no in Vc vector \n"
    << "        parametersVc[parameterNo*nVcVectors + i][k] = parameters[std::min(parameterNo*nInstances + i*(int)Vc::double_v::Size+k, nParametersPerInstance*nInstances-1)];\n"
    << std::endl
    << "  for (int i = 0; i < nVcVectors; i++)" << std::endl
    << "  {" << std::endl;

  // loop over lines of cellml code
  for (code_expression_t &codeExpression : cellMLCode_.lines)
  {
    if (codeExpression.type != code_expression_t::commented_out)
    {
      sourceCode << "    ";

      codeExpression.visitLeafs([&sourceCode,&nVcVectors,this](CellmlSourceCodeGeneratorVc::code_expression_t &expression, bool isFirstVariable)
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
              if (expression.code == "states")
              {
                sourceCode << "statesVc[" << expression.arrayIndex * nVcVectors << "+i]";
              }
              else if (expression.code == "rates")
              {
                sourceCode << "ratesVc[" << expression.arrayIndex * nVcVectors << "+i]";
              }
              else if (expression.code == "algebraics")
              {
                sourceCode << "algebraicsVc[" << expression.arrayIndex * nVcVectors << "+i]";
              }
              else if (expression.code == "parameters")
              {
                sourceCode << "parametersVc[" << expression.arrayIndex * nVcVectors << "+i]";
              }
              else
              {
                LOG(FATAL) << "unhandled variable type \"" << expression.code << "\".";
              }
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
        }
      });
      sourceCode << std::endl;
    }
  }
  sourceCode << std::endl;

  sourceCode << "  }" << std::endl << std::endl
    << "  // store computed values back to pointers\n"
    << "  for (int rateNo = 0; rateNo < nStates; rateNo++)\n"
    << "    for (int i = 0; i < nVcVectors; i++)  // Vc vector no\n"
    << "      for (int k = 0; k < Vc::double_v::Size; k++)  // entry no in Vc vector \n"
    << "      {\n"
    << "        if (rateNo*nInstances + i*Vc::double_v::Size+k >= nStates*nInstances)\n"
    << "          continue;\n"
    << "        rates[rateNo*nInstances + i*Vc::double_v::Size+k] = ratesVc[rateNo*nVcVectors + i][k];\n"
    << "      }\n"
    << "\n"
    << "  for (int algebraicNo = 0; algebraicNo < nAlgebraics; algebraicNo++)\n"
    << "    for (int i = 0; i < nVcVectors; i++)  // Vc vector no\n"
    << "      for (int k = 0; k < Vc::double_v::Size; k++)  // entry no in Vc vector \n"
    << "      {\n"
    << "        if (algebraicNo*nInstances + i*Vc::double_v::Size+k >= nAlgebraics*nInstances)\n"
    << "          continue;\n"
    << "        algebraics[algebraicNo*nInstances + i*Vc::double_v::Size+k] = algebraicsVc[algebraicNo*nVcVectors + i][k];\n"
    << "      }\n"
    << "\n";


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

  std::stringstream s;
  s << "-lVc -I\"" << OPENDIHU_HOME << "/dependencies/vc/install/include\" "
    << "-L\"" << OPENDIHU_HOME << "/dependencies/vc/install/lib\" ";
  if (std::string(CXX_COMPILER_COMMAND) == std::string("g++"))
    s << "-std=c++14";
  additionalCompileFlags_ = s.str();
  compilerCommand_ = CXX_COMPILER_COMMAND;
  sourceFileSuffix_ = ".cpp";
}

void CellmlSourceCodeGeneratorVc::
generateSourceFileVcFastMonodomain(std::string outputFilename, bool approximateExponentialFunction)
{
  std::set<std::string> helperFunctions;   //< functions found in the CellML code that need to be provided, usually the pow2, pow3, etc. helper functions for pow(..., 2), pow(...,3) etc.

  // replace pow and ?: functions
  preprocessCode(helperFunctions);

  std::stringstream sourceCode;
  sourceCode << "#include <math.h>" << std::endl
    << "#include <Vc/Vc>" << std::endl
    << "#include <iostream> " << std::endl
    << cellMLCode_.header << std::endl
    << "using Vc::double_v; " << std::endl;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  sourceCode << std::endl << "/* This file was created by opendihu at " << StringUtility::timeToString(&tm)  //std::put_time(&tm, "%d/%m/%Y %H:%M:%S")
    << ".\n * It is designed for the FastMonodomainSolver.\n "
    << " */\n";

  VLOG(1) << "call defineHelperFunctions with helperFunctions: " << helperFunctions;

  // define helper functions
  sourceCode << defineHelperFunctions(helperFunctions, approximateExponentialFunction);

  // define initializeStates function
  sourceCode
    << "// set initial values for all states\n"
    << "#ifdef __cplusplus\n" << "extern \"C\"\n" << "#endif\n" << std::endl
    << "void initializeStates(Vc::double_v states[]) \n"
    << "{\n";

  for (int stateNo = 0; stateNo < this->nStates_; stateNo++)
  {
    sourceCode << "  states[" << stateNo << "] = " << statesInitialValues_[stateNo] << ";\n";
  }
  sourceCode << "}\n\n";

  // define compute0D which computes one Heun step
  sourceCode
    << "// compute one Heun step\n"
    << "#ifdef __cplusplus\n" << "extern \"C\"\n" << "#endif\n" << std::endl
    << "void compute0DInstance(Vc::double_v states[], std::vector<Vc::double_v> &parameters, double currentTime, double timeStepWidth, bool stimulate,\n"
    << "                       bool storeAlgebraicsForTransfer, std::vector<Vc::double_v> &algebraicsForTransfer, const std::vector<int> &algebraicsForTransferIndices, double valueForStimulatedPoint) \n"
    << "{\n"
    << "  // assert that Vc::double_v::Size is the same as in opendihu, otherwise there will be problems\n"
    << "  if (Vc::double_v::Size != " << Vc::double_v::Size << ")\n"
    << "  {\n"
    << "    std::cout << \"Fatal error in compiled library of source file \\\"" << outputFilename << "\\\", size of SIMD register in "
    << "compiled code (\" << Vc::double_v::Size << \") does not match opendihu code (" << Vc::double_v::Size << ").\" << std::endl;\n"
    << "    std::cout << \"Delete library such that it will be regenerated with the correct compile options!\" << std::endl;\n"
    << "    exit(1);\n"
    << "  }\n\n"
    << "  // define constants\n";
    
/*    << R"(  std::cout << "currentTime=" << currentTime << ", timeStepWidth=" << timeStepWidth << ", stimulate=" << stimulate << std::endl;)" << "\n" */
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

    sourceCode << "  const double " << constantAssignmentsLine << std::endl;
  }
  sourceCode << "\n"
    << "  // compute new rates, rhs(y_n)\n";

  // loop over lines of cellml code
  for (code_expression_t &codeExpression : cellMLCode_.lines)
  {
    if (codeExpression.type != code_expression_t::commented_out)
    {
      std::stringstream sourceCodeLine;
      bool isCommentedOut = false;
      
      codeExpression.visitLeafs([&sourceCodeLine,&isCommentedOut,this](CellmlSourceCodeGeneratorVc::code_expression_t &expression, bool isFirstVariable)
      {
        switch(expression.type)
        {
          case code_expression_t::variableName:

            if (expression.code == "CONSTANTS")
            {
              // constants only exist once for all instances
              sourceCodeLine << "constant" << expression.arrayIndex;
            }
            else
            {
              // all other variables (states, rates, algebraics, parameters) exist for every instance
              if (expression.code == "states")
              {
                sourceCodeLine << "states[" << expression.arrayIndex << "]";
              }
              else if (expression.code == "rates")
              {
                sourceCodeLine << "rate" << expression.arrayIndex;
              }
              else if (expression.code == "algebraics")
              {
                sourceCodeLine << "algebraic" << expression.arrayIndex;
              }
              else if (expression.code == "parameters")
              {
                sourceCodeLine << "parameters[" << expression.arrayIndex << "]";
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
            sourceCodeLine << "  // (not assigning to a parameter) " << expression.code;
            isCommentedOut = true;
            break;

          default:
            break;
        }
      });
      
      if (isCommentedOut)
      {
        sourceCode << "  " << sourceCodeLine.str() << std::endl;
      }
      else
      {
        sourceCode << "  const double_v " << sourceCodeLine.str() << std::endl;
      }
    }
  }
  sourceCode << "\n"
    << "  // algebraic step\n"
    << "  // compute y* = y_n + dt*rhs(y_n), y_n = state, rhs(y_n) = rate, y* = algebraicState\n";

  for (int stateNo = 0; stateNo < this->nStates_; stateNo++)
  {
    sourceCode << "  ";
    if (stateNo != 0)
      sourceCode << "const ";

    sourceCode << "double_v algebraicState" << stateNo << " = states[" << stateNo << "] + timeStepWidth*rate" << stateNo << ";\n";
  }
  sourceCode << "\n\n"
    << R"(
  // if stimulation, set value of Vm (state0)
  if (stimulate)
  {
    for (int i = 0; i < std::min(3,(int)Vc::double_v::Size); i++)
    {
      algebraicState0[i] = valueForStimulatedPoint;
    }
  }

  // compute new rates, rhs(y*)
)";

  // loop over lines of cellml code
  for (code_expression_t &codeExpression : cellMLCode_.lines)
  {
    if (codeExpression.type != code_expression_t::commented_out)
    {
      std::stringstream sourceCodeLine;
      bool isCommentedOut = false;
      
      codeExpression.visitLeafs([&sourceCodeLine,&isCommentedOut,this](CellmlSourceCodeGeneratorVc::code_expression_t &expression, bool isFirstVariable)
      {
        switch(expression.type)
        {
          case code_expression_t::variableName:

            if (expression.code == "CONSTANTS")
            {
              // constants only exist once for all instances
              sourceCodeLine << "constant" << expression.arrayIndex;
            }
            else
            {
              // all other variables (states, rates, algebraics, parameters) exist for every instance
              if (expression.code == "states")
              {
                sourceCodeLine << "algebraicState" << expression.arrayIndex;
              }
              else if (expression.code == "rates")
              {
                sourceCodeLine << "algebraicRate" << expression.arrayIndex;
              }
              else if (expression.code == "algebraics")
              {
                sourceCodeLine << "algebraicAlgebraic" << expression.arrayIndex;
              }
              else if (expression.code == "parameters")
              {
                sourceCodeLine << "parameters[" << expression.arrayIndex << "]";
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
            sourceCodeLine << "  // (not assigning to a parameter) " << expression.code;
            isCommentedOut = true;
            break;

          default:
            break;
        }
      });
      
      if (isCommentedOut)
      {
        sourceCode << "  " << sourceCodeLine.str() << std::endl;
      }
      else
      {
        sourceCode << "  const double_v " << sourceCodeLine.str() << std::endl;
      }
    }
  }
  sourceCode << std::endl;

  sourceCode << R"(
  // final step
  // y_n+1 = y_n + 0.5*[rhs(y_n) + rhs(y*)]
)";

  for (int stateNo = 0; stateNo < this->nStates_; stateNo++)
  {
    sourceCode << "  states[" << stateNo << "] += 0.5*timeStepWidth*(rate" << stateNo << " + algebraicRate" << stateNo << ");\n";
  }

  sourceCode << R"(

  if (stimulate)
  {
    for (int i = 0; i < std::min(3,(int)Vc::double_v::Size); i++)
    {
      states[0][i] = valueForStimulatedPoint;
    }
  }

  // store algebraics for transfer
  if (storeAlgebraicsForTransfer)
  {
    for (int i = 0; i < algebraicsForTransferIndices.size(); i++)
    {
      const int algebraic = algebraicsForTransferIndices[i];

      switch (algebraic)
      {
)";

  // loop over algebraics and generate code to copy the updated algebraic values to the algebraics
  for (int algebraicNo = 0; algebraicNo < this->nAlgebraics_; algebraicNo++)
  {
    // only of the algebraic was computed and not replaced by a parameter
    if (std::find(this->parametersUsedAsAlgebraic_.begin(), this->parametersUsedAsAlgebraic_.end(), algebraicNo)
       != this->parametersUsedAsAlgebraic_.end())
    {
      sourceCode << "        // case " << algebraicNo << ": is a parameter\n";
    }
    else
    {
      sourceCode << "        case " << algebraicNo << ":\n"
        << "          algebraicsForTransfer[i] = algebraicAlgebraic" << algebraicNo << ";\n"
        << "          break;\n";
    }
  }

  sourceCode << R"(
      }
    }
  }
}
)";

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

  std::stringstream s;
  s << "-lVc -I\"" << OPENDIHU_HOME << "/dependencies/vc/install/include\" "
    << "-L\"" << OPENDIHU_HOME << "/dependencies/vc/install/lib\" ";
  if (std::string(CXX_COMPILER_COMMAND) == std::string("g++"))
    s << "-std=c++14";
  additionalCompileFlags_ = s.str();
  compilerCommand_ = CXX_COMPILER_COMMAND;
  sourceFileSuffix_ = ".cpp";
}

