#include "cellml/source_code_generator/source_code_generator.h"

#include <Python.h>  // has to be the first included header

#include "utility/string_utility.h"
#include "output_writer/generic.h"

#include <vector>
#include <iostream>
#include "easylogging++.h"

void CellMLSourceCodeGenerator::
parseSourceCodeFile()
{
  // input source filename is this->sourceFilename_

  // parse source file, set initial values for states (only one instance) and nParameters_, nConstants_ and nIntermediatesInSource_
  nIntermediatesInSource_ = 0;
  unsigned int nStatesInSource = 0;
  bool errorWrongNumberOfIntermediatesOrStates = false;
  bool currentlyInInitConstsFunction = false;

  // read in source from file
  std::ifstream sourceFile(this->sourceFilename_.c_str());
  if (!sourceFile.is_open())
  {
    LOG(FATAL) << "Could not open source file \"" << this->sourceFilename_<< "\" for reading!";
  }
  else
  {
    // read whole file contents to source
    std::stringstream source;
    source << sourceFile.rdbuf();
    sourceFile.close();

    bool discardOpenBrace = false;   // if the next line consisting of only "{" should be discarded
    bool headerDone = false;         // if everything before computeRates was parsed

    std::string name;  // the parsed name of a specifier that follows

    // step through lines and create simd source file
    while(!source.eof())
    {
      std::string line;
      getline(source, line);

      // discard the line with only "{"
      if (discardOpenBrace && line == "{")
      {
        discardOpenBrace = false;
        continue;
      }

      // put "void" on beginning of next line
      if (line == "void")
      {
        std::string line2;
        getline(source, line2);
        line += std::string(" ") + line2;
      }

      // parse initial values for states, constants and intermediates
      if (!headerDone)
      {
        if (line.find(" * STATES") == 0)  // line in OpenCOR generated input file of type " * STATES[55] is P_C_SR in component razumova (milliM)."
        {
          // parse name of state
          unsigned int index = atoi(line.substr(10,line.find("]")-10).c_str());
          int posBegin = line.find("is",12)+3;
          int posEnd = line.rfind(" in");
          name = line.substr(posBegin,posEnd-posBegin);
          VLOG(1) << "index= " << index << ", this->stateNames_.size() = " << this->stateNames_.size();
          nStatesInSource = std::max(nStatesInSource, index+1);

          if (index >= this->stateNames_.size())
          {
            LOG(ERROR) << "The CellML file \"" << sourceFilename_ << "\" contains more than " << index << " states "
              << " but only " << this->stateNames_.size() << " were given as template argument to CellMLAdapter.";
            errorWrongNumberOfIntermediatesOrStates = true;
          }
          else
          {
            this->stateNames_[index] = name;
          }
        }
        else if (line.find(" * ALGEBRAIC") == 0)  // line in OpenCOR generated input file of type " * ALGEBRAIC[35] is g_Cl in component sarco_Cl_channel (milliS_per_cm2)."
        {
          // parse name of intermediate
          unsigned int index = atoi(line.substr(13,line.find("]")-13).c_str());
          int posBegin = line.find("is",15)+3;
          int posEnd = line.rfind(" in");
          name = line.substr(posBegin,posEnd-posBegin);
          VLOG(1) << "index= " << index << ", this->intermediateNames_.size() = " << this->intermediateNames_.size();
          nIntermediatesInSource_ = std::max(nIntermediatesInSource_, index+1);

          if (index >= this->intermediateNames_.size())
          {
            LOG(ERROR) << "The CellML file \"" << sourceFilename_ << "\" contains more than " << index << " intermediates "
              << " but only " << this->intermediateNames_.size() << " were given as template argument to CellMLAdapter.";
            errorWrongNumberOfIntermediatesOrStates = true;
          }
          else
          {
            this->intermediateNames_[index] = name;
          }
        }
        else if (line.find("STATES[") == 0)   // line contains assignment in OpenCOR generated input file
        {
          // parse initial value of state
          unsigned int index = atoi(line.substr(7,line.find("]",7)-7).c_str());
          double value = atof(line.substr(line.find("= ")+2).c_str());
          if (index >= statesInitialValues_.size())
          {
            errorWrongNumberOfIntermediatesOrStates = true;
          }
          else
          {
            statesInitialValues_[index] = value;
          }
          nStatesInSource = std::max(nStatesInSource, index+1);
        }
        else if (line.find("ALGEBRAIC[") == 0)  // assignment to an algebraic variable in both OpenCMISS and OpenCOR generated files, in OpenCMISS generated files, this does not count towards the algebraic variables that are hold by opendihu
        {
          unsigned int algebraicIndex = atoi(line.substr(10,line.find("]",10)-10).c_str());
          nIntermediatesInSource_ = std::max(nIntermediatesInSource_, algebraicIndex+1);
        }
        else if (line.find("CONSTANTS[") != std::string::npos)  // usage of a constant
        {
          std::string substr(line.substr(line.find("CONSTANTS[")+10,line.find("]",line.find("CONSTANTS[")+10)-line.find("CONSTANTS[")-10));
          unsigned int index = atoi(substr.c_str());
          this->nConstants_ = std::max(this->nConstants_, index+1);
        }
        else if (line.find("OC_CellML_RHS_routine") != std::string::npos)
        {
          LOG(FATAL) << "Cellml sourceFilename \"" << sourceFilename_ << "\" is OpenCMISS generated. This is not supported. "
            << "Please use OpenCOR to generate the c source file.";
        }
      }

      // line contains declaration of ALGEBRAIC variable, e.g. "double CONSTANTS[110], ALGEBRAIC[70];"
      if (line.find("initConsts") != std::string::npos)
      {
        currentlyInInitConstsFunction = true;
      }
      else if (line.find("}") != std::string::npos && currentlyInInitConstsFunction)
      {
        currentlyInInitConstsFunction = false;
      }
      else if (line.find("CONSTANTS[") == 0)
      {
        constantAssignments_.push_back(line);
      }
      else if (line.find("void computeVariables") != std::string::npos)
      {
        // this is the end of the normal c source file
        break;
      }
      // line contains OpenCOR function head (computeRates)
      else if (line.find("computeRates") != std::string::npos)
      {
        discardOpenBrace = true;
        headerDone = true;
      }
      // line contains normal assignment
      else if (line.find("ALGEBRAIC") == 0 || line.find("RATES") == 0)
      {
        // parse line
        cellMLCode_.lines.emplace_back();

        cellMLCode_.lines.back().parse(line);

        VLOG(1) << "parse line [" << line << "] -> " << cellMLCode_.lines.back().getString();

        cellMLCode_.lines.back().visitLeafs([this](code_expression_t &entry, bool isFirstVariable)
        {
          if (entry.type == code_expression_t::variableName)
          {

            // check if this is an assignment to a algebraic value that is actually an explicit parameter (set by parametersUsedAsIntermediate)
            if (entry.code == "intermediates" && isFirstVariable)
            {
              for (int parameterUsedAsIntermediate : this->parametersUsedAsIntermediate_)
              {
                if (entry.arrayIndex == parameterUsedAsIntermediate)
                {
                  entry.type = code_expression_t::commented_out;
                  break;
                }
              }
            }

            // replace algebraic by parameter if it is an explicit parameter set by parametersUsedAsIntermediate
            else if (entry.code == "intermediates")
            {
              // loop over all parametersUsedAsIntermediate_
              for (int j = 0; j < this->parametersUsedAsIntermediate_.size(); j++)
              {
                if (entry.arrayIndex == this->parametersUsedAsIntermediate_[j])
                {
                  entry.code = "parameters";
                  entry.arrayIndex = j;
                  break;
                }
              }
            }

            // replace constant by parameter if it is an explicit parameter set by parametersUsedAsConstant
            else if (entry.code == "CONSTANTS")
            {
              // loop over all parametersUsedAsConstant_
              for (int j = 0; j < this->parametersUsedAsConstant_.size(); j++)
              {
                if (entry.arrayIndex == this->parametersUsedAsConstant_[j])
                {
                  entry.code = "parameters";
                  entry.arrayIndex = this->parametersUsedAsIntermediate_.size() + j;
                  break;
                }
              }
            }
          }
        });
      }
      // every other line
      else
      {
        if (!headerDone)
        {
          if (!currentlyInInitConstsFunction)
          {
            cellMLCode_.header += line + std::string("\n");
          }
        }
        else
          cellMLCode_.footer += line + std::string("\n");

        VLOG(2) << "line is not special, copy: [" << line << "]";
      }
    }

    if (errorWrongNumberOfIntermediatesOrStates)
    {
      LOG(FATAL) << "The CellML model \"" << sourceFilename_ << "\" has " << nStatesInSource << " states and "
        << nIntermediatesInSource_ << " intermediates, but the CellMLAdapter only supports "
        << nStates_ << " states and " << nIntermediates_ << " intermediates." << std::endl
        << "You have to set the correct number in the c++ file and recompile. " << std::endl
        << "(Use \" CellmlAdapter<" << nStatesInSource << "," << nIntermediatesInSource_ << ">\".)";
    }

    // check number of intermediates and states in source file
    if (nIntermediatesInSource_ != nIntermediates_)
    {
      LOG(WARNING) << "The CellML model \"" << sourceFilename_ << "\" needs " << nIntermediatesInSource_ << " intermediates and CellMLAdapter supports " << nIntermediates_
        << ". You should recompile with the correct number to avoid performance penalties." << std::endl
        << "(Use \" CellmlAdapter<" << nStatesInSource << "," << nIntermediatesInSource_ << ">\".)";
    }
    if (nStatesInSource != nStates_)
    {
      LOG(ERROR) << "The CellML model \"" << sourceFilename_ << "\" has " << nStatesInSource << " states and CellMLAdapter supports " << nStates_
        << ". This means the last " << nStates_ - nStatesInSource << " state(s) will have undefined values. You should recompile with the correct number of states." << std::endl
        << "(Use \" CellmlAdapter<" << nStatesInSource << "," << nIntermediatesInSource_ << ">\".)";
    }
  }
}

void CellMLSourceCodeGenerator::code_expression_t::
parse(std::string line)
{
  VLOG(2) << "line: [" << line << "]";

  std::vector<code_expression_t> entries;

  size_t currentPos = 0;
  for (int i = 0; currentPos < line.length(); i++)
  {
    VLOG(2);
    VLOG(2) << "currentPos: " << currentPos << " code from there: \"" << line.substr(currentPos, 20) << "\"";
    VLOG(2) << "variables (high number is string::npos and means not found): "
      << line.find("ALGEBRAIC", currentPos) << ", "
      << line.find("STATES", currentPos) << ", "
      << line.find("RATES", currentPos) << ", "
      << line.find("CONSTANTS", currentPos)
      << ", (=" << line.find("(", currentPos)
      << ", )=" << line.find(")", currentPos)
      << ", ?=" << line.find("?", currentPos);

    size_t posVariable = std::min({
      line.find("ALGEBRAIC", currentPos),   // algebraics
      line.find("STATES", currentPos),      // states
      line.find("RATES", currentPos),       // rates
      line.find("CONSTANTS", currentPos)    // CONSTANTS
    });

    std::size_t posOpeningParantheses = line.find("(", currentPos);
    std::size_t posClosingParantheses = line.find(")", currentPos);
    std::size_t posQuestionMark = line.find("?", currentPos);

    std::size_t posNextItem = std::min({
      posVariable,
      posOpeningParantheses,
      posClosingParantheses,
      posQuestionMark
    });

    if (currentPos < posNextItem)
    {
      code_expression_t entry;
      entry.type = code_expression_t::otherCode;
      entry.code = line.substr(currentPos, posNextItem - currentPos);
      entries.push_back(entry);

      currentPos = posNextItem;
      VLOG(2) << "parsed otherCode: [" << entry.code << "]";
    }
    else if (posOpeningParantheses == currentPos)
    {
      VLOG(2) << "opening parantheses, posClosingParantheses=" << posClosingParantheses;

      if (posClosingParantheses == std::string::npos)
      {
        LOG(FATAL) << "Missing closing paranthesis matching the opening paranthesis at pos " << posOpeningParantheses
          << " in line [" << line << "]";
      }

      // find matching closing parantheses
      int nOpenParantheses = 0;
      std::size_t i = posOpeningParantheses;
      for(;;)
      {
        std::size_t posNextOpeningParantheses = line.find("(", i);
        std::size_t posNextClosingParantheses = line.find(")", i);

        if (posNextOpeningParantheses < posNextClosingParantheses)
        {
          nOpenParantheses++;
          i = posNextOpeningParantheses + 1;
        }
        else
        {
          nOpenParantheses--;
          i = posNextClosingParantheses + 1;

          if (nOpenParantheses == 0)
          {
            posClosingParantheses = posNextClosingParantheses;
            break;
          }
        }
      }

      code_expression_t entry;
      entry.type = code_expression_t::tree;

      code_expression_t openingParanthesis;
      openingParanthesis.type = code_expression_t::otherCode;
      openingParanthesis.code = "(";

      code_expression_t body;
      std::string bodyString = line.substr(posOpeningParantheses+1, posClosingParantheses-(posOpeningParantheses+1));
      VLOG(2) << "bodyString: [" << bodyString << "]";

      body.parse(bodyString);

      code_expression_t closingParanthesis;
      closingParanthesis.type = code_expression_t::otherCode;
      closingParanthesis.code = ")";

      entry.treeChildren.push_back(openingParanthesis);
      entry.treeChildren.push_back(body);
      entry.treeChildren.push_back(closingParanthesis);
      entries.push_back(entry);

      currentPos = posClosingParantheses+1;
    }
    else if (posClosingParantheses == currentPos)
    {
      LOG(FATAL) << "Unmatched closing paranthesis at pos " << posClosingParantheses
          << " in line [" << line << "]";
    }
    else if (posQuestionMark == currentPos)
    {
      std::size_t posColon = line.find(":", currentPos);

      if (posColon == std::string::npos)
      {
        LOG(FATAL) << "Question mark does not have matching colon, at pos " << posQuestionMark
          << " in line [" << line << "]";
      }

      // find matching closing colon
      int nOpenTernaryOperators = 0;
      std::size_t i = posQuestionMark;
      for(;;)
      {
        std::size_t posNextQuestionMark = line.find("?", i);
        std::size_t posNextColon = line.find(":", i);

        if (posNextQuestionMark < posNextColon)
        {
          nOpenTernaryOperators++;
          i = posNextQuestionMark + 1;
        }
        else
        {
          nOpenTernaryOperators--;
          i = posNextColon + 1;

          if (nOpenTernaryOperators == 0)
          {
            posColon = posNextColon;
            break;
          }
        }
      }


      code_expression_t entry;
      entry.type = code_expression_t::tree;

      // condition of ternary operator
      code_expression_t condition;
      if (entries.size() == 1)
      {
        condition = entries[0];
      }
      else
      {
        condition.type = code_expression_t::tree;
        condition.treeChildren.assign(entries.begin(), entries.end());
        entries.clear();
      }

      // question mark
      code_expression_t questionMark;
      questionMark.type = code_expression_t::otherCode;
      questionMark.code = "?";

      // first branch
      code_expression_t firstBranch;
      std::string firstBranchString = line.substr(posQuestionMark+1, posColon-(posQuestionMark+1));
      VLOG(2) << "firstBranchString: [" << firstBranchString << "]";

      firstBranch.parse(firstBranchString);

      // colon
      code_expression_t colon;
      colon.type = code_expression_t::otherCode;
      colon.code = ":";

      // second branch
      code_expression_t secondBranch;
      std::string secondBranchString = line.substr(posColon+1);
      VLOG(2) << "secondBranchString: [" << secondBranchString << "]";

      secondBranch.parse(secondBranchString);

      entry.treeChildren.push_back(condition);
      entry.treeChildren.push_back(questionMark);
      entry.treeChildren.push_back(firstBranch);
      entry.treeChildren.push_back(colon);
      entry.treeChildren.push_back(secondBranch);
      entries.push_back(entry);

      currentPos = line.length();
    }
    else if (posVariable == currentPos)
    {
      code_expression_t entry;
      entry.type = code_expression_t::variableName;

      // parse code (variable name)
      size_t posBracket = line.find("[", currentPos);
      entry.code = line.substr(currentPos, posBracket-currentPos);

      // rename variable
      if (entry.code == "STATES")
        entry.code = "states";
      else if (entry.code == "RATES")
        entry.code = "rates";
      else if (entry.code == "ALGEBRAIC")
        entry.code = "intermediates";

      // extract array index
      entry.arrayIndex = atoi(line.substr(posBracket+1).c_str());

      VLOG(2) << "extract variable name \"" << entry.code << "\", index " << entry.arrayIndex;

      entries.push_back(entry);

      // advance current position
      currentPos = line.find("]", currentPos)+1;
    }
  }

  if (entries.size() == 1)
  {
    type = entries[0].type;
    code = entries[0].code;
    arrayIndex = entries[0].arrayIndex;
    treeChildren = entries[0].treeChildren;
  }
  else
  {
    type = code_expression_t::tree;
    treeChildren = entries;
  }
}


void CellMLSourceCodeGenerator::code_expression_t::
visitLeafs(std::function<void(CellMLSourceCodeGenerator::code_expression_t &expression, bool isFirstVariable)> callback)
{
  bool isFirstVariable = true;
  visitLeafsCall(callback, isFirstVariable);
}

void CellMLSourceCodeGenerator::code_expression_t::
visitLeafsCall(std::function<void(CellMLSourceCodeGenerator::code_expression_t &expression, bool isFirstVariable)> callback, bool &isFirstVariable)
{
  if (type == code_expression_t::tree)
  {
    // recursively iterate over children
    for (code_expression_t &innerExpression : treeChildren)
    {
      innerExpression.visitLeafsCall(callback, isFirstVariable);
    }
  }
  else
  {
    // call callback for own node
    callback(*this, isFirstVariable);

    // if this was a variable, set "isFirstVariable" to false
    if (type == code_expression_t::variableName)
      isFirstVariable = false;
  }
}

void CellMLSourceCodeGenerator::code_expression_t::
visitNodes(std::function<void(CellMLSourceCodeGenerator::code_expression_t &expression)> callback)
{
  // call callback for own node
  callback(*this);

  if (type == code_expression_t::tree)
  {
    // recursively iterate over children
    for (code_expression_t &innerExpression : treeChildren)
    {
      innerExpression.visitNodes(callback);
    }
  }
}

std::string CellMLSourceCodeGenerator::code_expression_t::getString()
{
  std::stringstream s;

  if (type == code_expression_t::tree)
  {
    s << "{";
    // recursively iterate over children
    for (code_expression_t &innerExpression : treeChildren)
    {
      s << innerExpression.getString() << ",";
    }
    s << "}";
  }
  else if (type == code_expression_t::variableName)
  {
    s << "'" << code << "[" << arrayIndex << "]'";
  }
  else if (type == code_expression_t::otherCode)
  {
    s << "'" << code << "'";
  }
  else if (type == code_expression_t::commented_out)
  {
    s << "'" << code << "'";
  }

  return s.str();
}
