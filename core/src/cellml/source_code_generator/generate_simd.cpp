#include "cellml/source_code_generator/source_code_generator.h"

#include <Python.h>  // has to be the first included header

#include "utility/string_utility.h"
#include "output_writer/generic.h"

#include <vector>
#include <iostream>
#include "easylogging++.h"

void CellMLSourceCodeGenerator::
generateSourceFileSimdAutovectorization(std::string outputFilename)
{
  std::stringstream simdSource;
  simdSource << "#include <math.h>" << std::endl
    << cellMLCode_.header << std::endl;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  simdSource << std::endl << "/* This function was created by opendihu at " << StringUtility::timeToString(&tm)  //std::put_time(&tm, "%d/%m/%Y %H:%M:%S")
    << ".\n * It is designed for " << this->nInstances_ << " instances of the CellML problem.\n "
    << "The \"optimizationType\" is \"simd\". (Other options are \"vc\" and \"openmp\".) */" << std::endl
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
  simdSource << std::endl
    << "  /* double ALGEBRAIC[" << this->nIntermediatesInSource_*this->nInstances_ << "];  */"
    << "        /* " << this->nIntermediatesInSource_ << " per instance * " << this->nInstances_ << " instances */ " << std::endl
    << "  double *ALGEBRAIC = intermediates;   /* use the storage for intermediates allocated by opendihu */" << std::endl << std::endl;

  // add all computations
  for (const std::list<assignment_t> &assignmentsList : cellMLCode_.entries)
  {
    if (assignmentsList.front().type != assignment_t::commented_out)
    {
      simdSource << std::endl
        << "#ifndef TEST_WITHOUT_PRAGMAS" << std::endl
        << "  #pragma omp for simd" << std::endl
        << "#endif" << std::endl
        //<< "#pragma GCC ivdep  // this disables alias checking for the compiler (GCC only)" << std::endl   // does not work on hazelhen
        << "  for (int i = 0; i < " << this->nInstances_ << "; i++)" << std::endl
        << "  {" << std::endl
        << "    ";

      VLOG(2) << "parsed " << cellMLCode_.entries.size() << " entries";

      // write out parsed code with adjusted indices
      for (auto entry : assignmentsList)
      {
        VLOG(2) << " entry type=" << entry.type << ", code: \"" << entry.code << "\", index: " << entry.arrayIndex;
        switch(entry.type)
        {
          case assignment_t::variableName:

            if (entry.code == "CONSTANTS")
            {
              // constants only exist once for all instances
              simdSource << entry.code << "[" << entry.arrayIndex<< "]";

              VLOG(2) << "    (write \"" << entry.code << "[" << entry.arrayIndex<< "]" << "\")";
            }
            else
            {
              // all other variables (states, rates, intermediates, parameters) exist for every instance
              simdSource << entry.code << "[" << entry.arrayIndex * this->nInstances_ << "+i]";

              VLOG(2) << "    (write \"" << entry.code << "[" << entry.arrayIndex * this->nInstances_ << "+i]" << "\")";
            }
            break;

          case assignment_t::other:

            VLOG(2) << "    (write \"" << entry.code << "\")";

            simdSource << entry.code;
            break;

          case assignment_t::commented_out:

            simdSource << "  // (not assigning to a parameter) " << entry.code;
            break;
        }
      }
      VLOG(2) << "write end of for loop (closing })";
      simdSource << std::endl << "  }" << std::endl;
    }
  }

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
}

void CellMLSourceCodeGenerator::
generateSourceFileExplicitVectorization(std::string outputFilename)
{
  LOG(WARNING) << "In CellmlAdapter optimizationType \"vc\" is not yet implemented, fallback to \"simd\".";
  generateSourceFileSimdAutovectorization(outputFilename);
}

void CellMLSourceCodeGenerator::
generateSourceFileOpenMP(std::string outputFilename)
{
  std::stringstream simdSource;
  simdSource << "#include <math.h>" << std::endl
    << "#include <omp.h>" << std::endl
    << cellMLCode_.header << std::endl;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  simdSource << std::endl << "/* This function was created by opendihu at " << StringUtility::timeToString(&tm)  //std::put_time(&tm, "%d/%m/%Y %H:%M:%S")
    << ".\n * It is designed for " << this->nInstances_ << " instances of the CellML problem.\n "
    << "The \"optimizationType\" is \"openmp\". (Other options are \"vc\" and \"simd\".) */" << std::endl
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
  simdSource << std::endl
    << "  /* double ALGEBRAIC[" << this->nIntermediatesInSource_*this->nInstances_ << "];  */"
    << "        /* " << this->nIntermediatesInSource_ << " per instance * " << this->nInstances_ << " instances */ " << std::endl
    << "  double *ALGEBRAIC = intermediates;   /* use the storage for intermediates allocated by opendihu */" << std::endl << std::endl;

  simdSource << std::endl
    << "  #pragma omp parallel for" << std::endl
    << "  for (int i = 0; i < " << this->nInstances_ << "; i++)" << std::endl
    << "  {" << std::endl;

  // add all computations
  for (const std::list<assignment_t> &assignmentsList : cellMLCode_.entries)
  {
    if (assignmentsList.front().type != assignment_t::commented_out)
    {
      VLOG(2) << "parsed " << cellMLCode_.entries.size() << " entries";

      simdSource << "    ";

      // write out parsed code with adjusted indices
      for (auto entry : assignmentsList)
      {
        VLOG(2) << " entry type=" << entry.type << ", code: \"" << entry.code << "\", index: " << entry.arrayIndex;
        switch(entry.type)
        {
          case assignment_t::variableName:

            if (entry.code == "CONSTANTS")
            {
              // constants only exist once for all instances
              simdSource << entry.code << "[" << entry.arrayIndex<< "]";

              VLOG(2) << "    (write \"" << entry.code << "[" << entry.arrayIndex<< "]" << "\")";
            }
            else
            {
              // all other variables (states, rates, intermediates, parameters) exist for every instance
              simdSource << entry.code << "[" << entry.arrayIndex * this->nInstances_ << "+i]";

              VLOG(2) << "    (write \"" << entry.code << "[" << entry.arrayIndex * this->nInstances_ << "+i]" << "\")";
            }
            break;

          case assignment_t::other:

            VLOG(2) << "    (write \"" << entry.code << "\")";

            simdSource << entry.code;
            break;

          case assignment_t::commented_out:

            simdSource << "  // (not assigning to a parameter) " << entry.code;
            break;
        }
      }

      simdSource << std::endl;
    }
  }
  VLOG(2) << "write end of for loop (closing })";
  simdSource << std::endl << "  }" << std::endl;

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

  additionalCompileFlags_ = "-fopenmp";
}

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
          LOG(DEBUG) << "index= " << index << ", this->stateNames_.size() = " << this->stateNames_.size();
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
          LOG(DEBUG) << "index= " << index << ", this->intermediateNames_.size() = " << this->intermediateNames_.size();
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
        cellMLCode_.entries.emplace_back();
        std::list<assignment_t> &entries = cellMLCode_.entries.back();

        bool isExplicitParameter = false;   // if this is an explicit parameter, i.e. a ALGEBRAIC variable, that is overridden by a parameter

        VLOG(2) << "line: [" << line << "]";

        size_t currentPos = 0;
        for (int i = 0; currentPos <= line.length(); i++)
        {
          VLOG(2);
          VLOG(2) << "currentPos: " << currentPos << " code from there: \"" << line.substr(currentPos, 20) << "\"";
          VLOG(2) << "variables (high number is string::npos and means not found): "
            << line.find("ALGEBRAIC", currentPos) << ", "
            << line.find("STATES", currentPos) << ", "
            << line.find("RATES", currentPos) << ", "
            << line.find("CONSTANTS", currentPos);

          size_t posVariable = std::min({
            line.find("ALGEBRAIC", currentPos),   // algebraics
            line.find("STATES", currentPos),      // states
            line.find("RATES", currentPos),       // rates
            line.find("CONSTANTS", currentPos)    // CONSTANTS
          });

          VLOG(2) << "posVariable: " <<posVariable;

          assignment_t entry;
          if (posVariable == currentPos)
          {
            size_t posBracket = line.find("[", currentPos);
            entry.code = line.substr(currentPos, posBracket-currentPos);

            // rename variable
            if (entry.code == "STATES")
              entry.code = "states";
            else if (entry.code == "RATES")
              entry.code = "rates";
            else if (entry.code == "ALGEBRAIC")
              entry.code = "ALGEBRAIC";
            
            // extract array index
            entry.arrayIndex = atoi(line.substr(posBracket+1).c_str());

            VLOG(2) << "extract variable name \"" << entry.code << "\", index " << entry.arrayIndex;

            entry.type = assignment_t::variableName;

            // check if number of states is in bounds of template parameter nStates
            if (entry.code == "states" && entry.arrayIndex >= nStates_)
            {
              LOG(FATAL) << "CellML code in source file \"" << this->sourceFilename_ << "\" "
                << "computes more states than given in the user code as template parameter "
                << "(template parameter: " << nStates_ << ", CellML code: at least " << entry.arrayIndex+1 << ").";
            }

            // check if this is an assignment to a algebraic value that is actually an explicit parameter (set by parametersUsedAsIntermediate)
            if (entry.code == "ALGEBRAIC" && i == 0)
            {
              isExplicitParameter = false;
              for (int parameterUsedAsIntermediate : this->parametersUsedAsIntermediate_)
              {
                if (entry.arrayIndex == parameterUsedAsIntermediate)
                {
                  isExplicitParameter = true;
                  break;
                }
              }

              if (isExplicitParameter)
              {
                entry.type = assignment_t::commented_out;
                entry.code = line;
                entries.push_back(entry);
                entries.front().type = assignment_t::commented_out;
                break;
              }
            }

            // replace algebraic by parameter if it is an explicit parameter set by parametersUsedAsIntermediate
            if (entry.code == "ALGEBRAIC")
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
            if (entry.code == "CONSTANTS")
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

            // advance current position
            currentPos = line.find("]", currentPos)+1;
            VLOG(2) << "(1)advance to " << currentPos;
          }
          else
          {
            entry.code = line.substr(currentPos, posVariable-currentPos);
            VLOG(2) << "extract code \"" << entry.code << "\".";
            entry.type = assignment_t::other;

            // advance current position
            currentPos = posVariable;
            VLOG(2) << "(2)advance to " << currentPos;
          }
          entries.push_back(entry);
        }

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
