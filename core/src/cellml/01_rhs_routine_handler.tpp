#include "cellml/01_rhs_routine_handler.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/mesh_manager.h"

#include <unistd.h>  //dlopen
#include <dlfcn.h> 
#include <ctime>

// forward declaration
template <int nStates>
class CellmlAdapter;


template<int nStates>
void RhsRoutineHandler<nStates>::
initializeRhsRoutine()
{
  // 1) if simdSourceFilename is given, use that source to compile the library
  // 2) if not 1) but sourceFilename is given, create simdSourceFilename from that and compile library
  // 3) if not 2) but libraryFilename is given, load that library, if it contains simdRhs, use that, if it contains non-simd rhs use that

  // try to load or create simd source
  std::string simdSourceFilename;
  if (PythonUtility::hasKey(this->specificSettings_, "sourceFilename"))
  {
    if (!createSimdSourceFile(simdSourceFilename))
      simdSourceFilename = "";
  }

  if (simdSourceFilename == "")
  {
    if(PythonUtility::hasKey(this->specificSettings_, "simdSourceFilename"))
    {
      simdSourceFilename = PythonUtility::getOptionString(this->specificSettings_, "simdSourceFilename", "");
    }
  }

  std::string libraryFilename;
  // if simdSourceFilename is set, compile to create dynamic library
  if (simdSourceFilename != "")
  {
    // compile source file to a library
    libraryFilename = "lib.so";
    if (PythonUtility::hasKey(this->specificSettings_, "libraryFilename"))
    {
      libraryFilename = PythonUtility::getOptionString(this->specificSettings_, "libraryFilename", "lib.so");
    }
    else 
    {
      std::stringstream s;
      s << "lib/"+StringUtility::extractBasename(this->sourceFilename_) << "_" << this->nInstances_ << ".so";
      libraryFilename = s.str();
    }

    bool doCompilation = true;
    forceRecompileRhs_ = PythonUtility::getOptionBool(this->specificSettings_, "forceRecompileRhs", true);
    if (!forceRecompileRhs_)
    {
      // check if the library file already exists
      std::ifstream file;
      file.open(libraryFilename);
      if (file.is_open())
      {
        LOG(DEBUG) << "Library \"" << libraryFilename << "\" already exists, do not recompile (set forceRecompileRhs to True to force recompilation).";
        doCompilation = false;
        file.close();
      }
    }
    
    if (doCompilation)
    {
     
      if (libraryFilename.find("/") != std::string::npos)
      {
        std::string path = libraryFilename.substr(0, libraryFilename.rfind("/"));
        int ret = system((std::string("mkdir -p ")+path).c_str());
        
        if (ret != 0)
        {
          LOG(ERROR) << "Could not create path \"" << path << "\".";
        }
      }
     
      std::stringstream compileCommand;
      // -ftree-vectorize -fopt-info-vec-missed -fopt-info-vec-optimized
#ifdef NDEBUG
      compileCommand << "gcc -fPIC -O3 -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared -lm -x c -o " << libraryFilename << " " << simdSourceFilename;
#else
      compileCommand << "gcc -fPIC -O0 -ggdb -shared -lm -x c -o " << libraryFilename << " " << simdSourceFilename;
#endif

      int ret = system(compileCommand.str().c_str());
      if (ret != 0)
      {
        LOG(ERROR) << "Compilation failed. Command: \"" << compileCommand.str() << "\".";
        libraryFilename = "";
      }
      else
      {
        LOG(DEBUG) << "Compilation successful. Command: \"" << compileCommand.str() << "\".";
      }

      // repeat compilation with different GCC vectorizer outputs
  #if 0
      compileCommand.str("");
      compileCommand << "gcc -fPIC -O3 -ftree-vectorize -fopt-info-vec-missed=vectorizer_missed.log -shared -lm -x c -o " << libraryFilename << " " << simdSourceFilename;
      ret = system(compileCommand.str().c_str());
      if (ret != 0)
      {
        LOG(DEBUG) << "Compilation failed. Command: \"" << compileCommand.str() << "\".";
      }
      else
      {
        LOG(DEBUG) << "Compilation successful. Command: \"" << compileCommand.str() << "\".";
      }
      compileCommand.str("");
      compileCommand << "gcc -fPIC -O3 -ftree-vectorize -fopt-info-vec-all=vectorizer_all.log -shared -lm -x c -o " << libraryFilename << " " << simdSourceFilename;
      ret = system(compileCommand.str().c_str());
      if (ret != 0)
      {
        LOG(DEBUG) << "Compilation failed. Command: \"" << compileCommand.str() << "\".";
      }
      else
      {
        LOG(DEBUG) << "Compilation successful. Command: \"" << compileCommand.str() << "\".";
      }
  #endif
    }
  }

  // if library is still not available, look if it was provided
  if (libraryFilename == "")
  {
    libraryFilename = PythonUtility::getOptionString(this->specificSettings_, "libraryFilename", "lib.so");
  }

  if (libraryFilename == "")
  {
    LOG(FATAL) << "Could not create or locate dynamic library for cellml right hand side routine.";
  }

  loadRhsLibrary(libraryFilename);
}

template<int nStates>
bool RhsRoutineHandler<nStates>::
loadRhsLibrary(std::string libraryFilename)
{
 
  // load dynamic library
  //
  std::string currentWorkingDirectory = getcwd(NULL,0);
  // append slash if there is none at the end
  if (currentWorkingDirectory[currentWorkingDirectory.length()-1] != '/')
    currentWorkingDirectory += "/";

  void* handle = dlopen((currentWorkingDirectory+libraryFilename).c_str(), RTLD_LOCAL | RTLD_LAZY);
  if (handle)
  {
    // load rhs method
    
    // try to load several routine names
    rhsRoutine_ = (void (*)(void *,double,double*,double*,double*,double*)) dlsym(handle, "computeCellMLRightHandSide");
    rhsRoutineOpenCMISS_ = (void (*)(double,double*,double*,double*,double*)) dlsym(handle, "OC_CellML_RHS_routine");
    initConstsOpenCOR_ = (void(*)(double*, double*, double*)) dlsym(handle, "initConsts");
    computeRatesOpenCOR_ = (void(*)(double, double*, double*, double*, double*)) dlsym(handle, "computeRates");
    computeVariablesOpenCOR_  = (void(*)(double, double*, double*, double*, double*)) dlsym(handle, "computeVariables");

    LOG(DEBUG) << "Library \"" << libraryFilename << "\" loaded. "
      << "rhsRoutine_: " << (rhsRoutine_==NULL? "NULL" : "yes")
      << ", rhsRoutineOpenCMISS_: " << (rhsRoutineOpenCMISS_==NULL? "NULL" : "yes")
      << ", initConstsOpenCOR_: " << (initConstsOpenCOR_==NULL? "NULL" : "yes")
      << ", rhsRoutineOpenCMISS_: " << (rhsRoutineOpenCMISS_==NULL? "NULL" : "yes")
      << ", computeRatesOpenCOR_: " << (computeRatesOpenCOR_==NULL? "NULL" : "yes")
      << ", computeVariablesOpenCOR_: " << (computeVariablesOpenCOR_==NULL? "NULL" : "yes");

    if (rhsRoutine_)
    {
      // if the opendihu-generated rhs function (with openmp pragmas) is present in the library, we can directly use it and we are done in this method.
      return true;
    }
    else if (rhsRoutineOpenCMISS_)
    {
      // if the OpenCMISS routine is present, setup the rhs function such that the OpenCMISS function is called for every instance of the problem
      rhsRoutine_ = [](void *context, double t, double *statesInput, double *ratesOutput, double *algebraicsOutput, double *parametersInput)
      {
        LOG(DEBUG) << "call opendihu rhsRoutine, by calling several opencmiss rhs";

        CellmlAdapter<nStates> *cellmlAdapter = (CellmlAdapter<nStates> *)context;
        int nInstances, nIntermediates, nParameters;
        cellmlAdapter->getNumbers(nInstances, nIntermediates, nParameters);

        // call the standard rhs routine for every instance of the CellML problem
        for(int instanceNo = 0; instanceNo < nInstances; instanceNo++)
        {
          // create vectors that contain the values for a single instance
          std::vector<double> states(nStates);
          std::vector<double> rates(nStates);
          std::vector<double> intermediates(nIntermediates);
          std::vector<double> parameters(nParameters);

          // copy the values from the interleaved vectors
          for (int i = 0; i < nStates; i++)
          {
            states[i] = statesInput[i*nInstances + instanceNo];
          }
          for (int i = 0; i < nParameters; i++)
          {
            parameters[i] = parametersInput[i*nInstances + instanceNo];
          }

          // call OpenCMISS generated rhs routine
          cellmlAdapter->rhsRoutineOpenCMISS_(t, states.data(), rates.data(), intermediates.data(), parameters.data());

          // copy the results back in the big vectors

          // copy the values back to the interleaved vectors
          for (int i = 0; i < nStates; i++)
          {
            ratesOutput[i*nInstances + instanceNo] = rates[i];
          }

          for (int i = 0; i < nIntermediates; i++)
          {
            algebraicsOutput[i*nInstances + instanceNo] = intermediates[i];
          }
        }
      };
      
    }
    else if (initConstsOpenCOR_)
    {
      LOG(FATAL) << "Could not load rhs routine from dynamic library \"" << libraryFilename << "\", but a computeRates method from OpenCOR was detected. " 
        << "Please provide the sourceFilename of the source instead of the library.";
      return false;
    }
    else
    {
      LOG(FATAL) << "Could not load rhs routine from dynamic library \"" << libraryFilename << "\".";
      return false;
    }
  }
  else
  {
    LOG(FATAL) << "Could not load dynamic library \"" << libraryFilename << "\". Reason: " <<dlerror();
  } 
}

template<int nStates>
bool RhsRoutineHandler<nStates>::
createSimdSourceFile(std::string &simdSourceFilename)
{
  // This method can handle two different types of input c files: from OpenCMISS and from OpenCOR
  
  // input source filename is this->sourceFilename_
  bool inputFileTypeOpenCMISS = true;   //< if the input file that is being parsed is from OpenCMISS and not from OpenCOR
  
  // read in source from file
  std::ifstream sourceFile(this->sourceFilename_.c_str());
  if (!sourceFile.is_open())
  {
    LOG(ERROR) << "Could not open source file \"" << this->sourceFilename_<< "\" for reading!";
    return false;
  }
  else
  {
    // read whole file contents to source
    std::stringstream source;
    source << sourceFile.rdbuf();
    sourceFile.close();

    bool discardOpenBrace = false;   // if the next line consisting of only "{" should be discarded
    std::stringstream simdSource;
    
    simdSource << "#include <math.h>" << std::endl;
    
    // step through lines and create simd source file
    while(!source.eof())
    {
      std::string line;
      getline(source, line);
      
      if (discardOpenBrace && line == "{")
      {
        discardOpenBrace = false;
        continue;
      }
      
      if (line == "void")
      {
        std::string line2;
        getline(source, line2);
        line += std::string(" ") + line2;
      }

      // line contains declaration of ALGEBRAIC variable, e.g. "double CONSTANTS[110], ALGEBRAIC[70];"
      if (line.find("double CONSTANTS") == 0 && line.find("ALGEBRAIC[") != std::string::npos)
      {
        size_t pos = line.find("ALGEBRAIC[");   
        size_t posStart = line.find("[", pos)+1;
        size_t posEnd = line.find("]", posStart);
        int algebraicSize = atoi(line.substr(posStart).c_str());

        simdSource << line.substr(0, posStart) << algebraicSize * this->nInstances_ << line.substr(posEnd) << std::endl;
      }
      else if (line.find("initConsts") != std::string::npos)
      {
        inputFileTypeOpenCMISS = false;
        simdSource << line << std::endl;
      }
      else if (line.find("CONSTANTS[") == 0)
      {
        constantAssignments_.push_back(line);
        simdSource << line << std::endl;
      }
      else if (line.find("void computeVariables") != std::string::npos)
      {
        simdSource << std::endl;
        break;
      }
      // line contains OpenCMISS function head
      else if (line.find("void OC_CellML_RHS_routine") != std::string::npos || line.find("computeRates") != std::string::npos)
      {
        if (line.find("void OC_CellML_RHS_routine") != std::string::npos)
        {
          inputFileTypeOpenCMISS = true;
        }
        else
        {
          inputFileTypeOpenCMISS = false; 
        }
        
        if (line.find("void computeCellMLRightHandSide") != std::string::npos)
        {
          LOG(WARNING) << "The given source file \"" << this->sourceFilename_<< "\" already contains a simd version of the rhs routine. "
            << "Use the option \"simdSourceFilename\" instead.";
          return true;
        }

        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        simdSource << std::endl << "/* This function was created by opendihu at " << std::put_time(&tm, "%d/%m/%Y %H:%M:%S") 
          << ".\n * It is designed for " << this->nInstances_ << " instances of the CellML problem. */" << std::endl
          << "void computeCellMLRightHandSide("
          << "void *context, double t, double *states, double *rates, double *algebraics, double *parameters)" << std::endl << "{" << std::endl;
        discardOpenBrace = true;
          
        simdSource << "  double VOI = t;   /* current simulation time */" << std::endl;
        if (!inputFileTypeOpenCMISS) 
        {
          simdSource << std::endl << "  /* define constants */" << std::endl 
            << "  double CONSTANTS[" << this->nConstants_ << "];" << std::endl;
          for (std::string constantAssignmentsLine : constantAssignments_)
          {
            simdSource << "  " << constantAssignmentsLine << std::endl; 
          }
          simdSource << std::endl
            << "  double ALGEBRAIC[" << this->nIntermediates_*this->nInstances_ << "];  " 
            << "  /* " << this->nIntermediates_ << " per instance * " << this->nInstances_ << " instances */ " << std::endl;
        }
      }
      // line contains OpenCMISS assignment
      else if(line.find("OC_WANTED") == 0 || line.find("OC_RATE") == 0 || line.find("ALGEBRAIC") == 0 || line.find("RATES") == 0)
      {
        // parse line
        struct entry_t
        {
          std::string code;
          int arrayIndex;
          enum {variableName, other} type;
        };
        std::list<entry_t> entries;

        bool isExplicitParameter;   // if this is an explicit parameter, i.e. a ALGEBRAIC variable, that is overridden by a parameter
          
        VLOG(2) << "line: [" << line << "]";

        size_t currentPos = 0;
        for(int i = 0; currentPos <= line.length(); i++)
        {
          VLOG(2);
          VLOG(2) << "currentPos: " << currentPos << " code from there: \"" << line.substr(currentPos, 20) << "\"";
          VLOG(2) << "variables: "
            << line.find("OC_STATE", currentPos) << ", "
            << line.find("OC_RATE", currentPos) << ", "
            << line.find("ALGEBRAIC", currentPos) << ", "
            << line.find("OC_WANTED", currentPos) << ", "
            << line.find("OC_KNOWN", currentPos) << ", "
            << line.find("STATES", currentPos) << ", "
            << line.find("RATES", currentPos) << ", "
            << line.find("CONSTANTS", currentPos);

          size_t posVariable = std::min({
            line.find("OC_STATE", currentPos),    // states
            line.find("OC_RATE", currentPos),     // rates
            line.find("ALGEBRAIC", currentPos),   // algebraics
            line.find("OC_WANTED", currentPos),   // algebraics
            line.find("OC_KNOWN", currentPos),    // parameters
            line.find("STATES", currentPos),      // states
            line.find("RATES", currentPos),       // rates
            line.find("CONSTANTS", currentPos)    // CONSTANTS
          });

          VLOG(2) << "posVariable: " <<posVariable;

          entry_t entry;
          if (posVariable == currentPos)
          {
            size_t posBracket = line.find("[", currentPos);
            entry.code = line.substr(currentPos, posBracket-currentPos);
            
            // rename variable
            if (entry.code == "OC_STATE" || entry.code == "STATES")
              entry.code = "states";
            else if (entry.code == "OC_RATE" || entry.code == "RATES")
              entry.code = "rates";
            else if (entry.code == "ALGEBRAIC" || entry.code == "OC_WANTED")
              entry.code = "algebraics";
            else if (entry.code == "OC_KNOWN")
              entry.code = "parameters";
            
            // extract array index
            entry.arrayIndex = atoi(line.substr(posBracket+1).c_str());

            VLOG(2) << "extract variable name \"" << entry.code << "\", index " << entry.arrayIndex;

            entry.type = entry_t::variableName;

            // check if number of states is in bounds of template parameter nStates
            if (entry.code == "states" && entry.arrayIndex >= nStates)
            {
              LOG(FATAL) << "CellML code in source file \"" << this->sourceFilename_ << "\" "
                << "computes more states than given in the user code as template parameter "
                << "(template parameter: " << nStates << ", CellML code: at least " << entry.arrayIndex+1 << ").";
            }
            
            // check if this is an assignment to a algebraic value that is actually an explicit parameter (set by parametersUsedAsIntermediate)
            if (entry.code == "algebraics" && i == 0)
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
                break;
            }
            
            // replace algebraic by parameter if it is an explicit parameter set by parametersUsedAsIntermediate
            if (entry.code == "algebraics")
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
            entry.type = entry_t::other;

            // advance current position
            currentPos = posVariable;
            VLOG(2) << "(2)advance to " << currentPos;
          }
          entries.push_back(entry);
        }
        
        if (isExplicitParameter)
        {
          simdSource << "  /* explicit parameter */" << std::endl 
            << "  /* " << line << "*/" << std::endl;
        }
        else
        { 
          // add pragma omp here
          simdSource << std::endl << "  for(int i = 0; i < " << this->nInstances_ << "; i++)" << std::endl
            << "  {" << std::endl << "    ";

          VLOG(2) << "parsed " << entries.size() << " entries";

          // write out parsed code with adjusted indices
          for (auto entry : entries)
          {
            VLOG(2) << " entry type=" << entry.type << ", code: \"" << entry.code << "\", index: " << entry.arrayIndex;
            switch(entry.type)
            {
              case entry_t::variableName:
                if (entry.code == "CONSTANTS")
                {
                  // constants only exist once for all instances
                  simdSource << entry.code << "[" << entry.arrayIndex<< "]";
                }
                else 
                {
                  // all other variables (states, rates, intermediates, parameters) exist for every instance
                  simdSource << entry.code << "[" << entry.arrayIndex * this->nInstances_ << "+i]";
                }
                break;
              case entry_t::other:
                simdSource << entry.code;
                break;
            }
          }
          simdSource << std::endl << "  }" << std::endl;
        }
      }
      // every other line
      else
      {
        simdSource << line << std::endl;
      }
    }

    // write out source file
    std::stringstream s;
    s << StringUtility::extractBasename(this->sourceFilename_) << "_simd.c";
    simdSourceFilename = s.str();
    if (PythonUtility::hasKey(this->specificSettings_, "simdSourceFilename"))
    {
      simdSourceFilename = PythonUtility::getOptionString(this->specificSettings_, "simdSourceFilename", "");
    }

    std::ofstream simdSourceFile(simdSourceFilename.c_str());
    if (!simdSourceFile.is_open())
    {
      LOG(ERROR) << "Could not write to file \"" << simdSourceFilename << "\".";
      return false;
    }
    else
    {
      std::string fileContents = simdSource.str();
      simdSourceFile << fileContents;
      simdSourceFile.close();
    }
  }
  
  return true;
}

template<int nStates>
bool RhsRoutineHandler<nStates>::
scanSourceFile(std::string sourceFilename, std::array<double,nStates> &statesInitialValues)
{
  LOG(TRACE) << "scanSourceFile";
  
  // parse source file, set initial values for states (only one instance) and nParameters_, nConstants_ and nIntermediates_

  this->inputFileTypeOpenCMISS_ = true;   //< if the input file that is being parsed is from OpenCMISS and not from OpenCOR
  
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

      // line contains initial value for a state or a known value, for example: "DUMMY_ASSIGNMENT /*OC_STATE[0]*/ = -79.974;", OpenCMISS
      if (line.find("DUMMY_ASSIGNMENT") == 0 && line.find("/*") != std::string::npos)
      {
        // parse line
        size_t posBegin = line.find("/*")+2;
        size_t posEnd = line.find("*/");
        std::string variableType = line.substr(posBegin, posEnd-posBegin);
        size_t pos = variableType.find("[");
        std::string variableName = variableType.substr(0, pos);
        unsigned int index = atoi(variableType.substr(pos+1).c_str());

        pos = line.find("= ");
        double value = atof(line.substr(pos+2).c_str());

        if (variableName == "OC_STATE" && index >= 0 && index < (unsigned int)nStates)
        {
          // store initial value for state
          statesInitialValues[index] = value;
          
          // store name of the state that was parsed from the comment in the previous line
          this->stateNames_[index] = name;
        }
        else if (variableName == "OC_KNOWN" && index >= 0)
        {
          // store initial value for parameter
          if (this->parameters_.size() < index+1)
          {
            this->parameters_.resize(index+1);
          }
          this->parameters_[index] = value;
        }
      }
      else if (line.find("/* Constant") == 0)  // line in OpenCMISS generated input file of type "/* Constant C_m */"
      {
        // parse name of state, will be assigned to the state in the next line
        name = line.substr(12, line.find(" */")-12);
      }
      else if (line.find(" * STATES") == 0)  // line in OpenCOR generated input file of type " * STATES[55] is P_C_SR in component razumova (milliM)."
      {
        // parse name of state
        unsigned int index = atoi(line.substr(10,line.find("]")-10).c_str());
        int posBegin = line.find("is",12)+3;
        int posEnd = line.rfind(" in");
        name = line.substr(posBegin,posEnd-posBegin);
        LOG(DEBUG) << "index= " << index << ", this->stateNames_size = " << this->stateNames_.size();
        this->stateNames_[index] = name;
      }
      else if (line.find("STATES[") == 0)   // line contains assignment in OpenCOR generated input file
      {
        // parse initial value of state
        unsigned int index = atoi(line.substr(7,line.find("]",7)-7).c_str());
        double value = atof(line.substr(line.find("= ")+2).c_str());
        statesInitialValues[index] = value;
      }
      else if (line.find("ALGEBRAIC[") == 0)  // assignment to an algebraic variable in both OpenCMISS and OpenCOR generated files
      {
        int algebraicIndex = atoi(line.substr(10,line.find("]",10)-10).c_str());
        this->nIntermediates_ = std::max(this->nIntermediates_, algebraicIndex+1);
      }
      else if (line.find("OC_KNOWN[") != std::string::npos)  // usage of a parameter variable in OpenCMISS generated file
      {
        std::string substr(line.substr(line.find("OC_KNOWN[")+9,line.find("]",line.find("OC_KNOWN[")+9)-line.find("OC_KNOWN[")-9));
        int index = atoi(substr.c_str());
        this->nParameters_ = std::max(this->nParameters_, index+1);
      }
      else if (line.find("CONSTANTS[") != std::string::npos)  // usage of a constant
      {
        std::string substr(line.substr(line.find("CONSTANTS[")+10,line.find("]",line.find("CONSTANTS[")+10)-line.find("CONSTANTS[")-10));
        int index = atoi(substr.c_str());
        this->nConstants_ = std::max(this->nConstants_, index+1);
      }
      
      if (line.find("initConsts") == 0)
      {
        this->inputFileTypeOpenCMISS_ = false;
      }
      
      // if the rhs routine is reached in the OpenCOR file, there are no more initializations, there quit processing of the source file
      if (!this->inputFileTypeOpenCMISS_ && line.find("computeCellMLRightHandSide") != std::string::npos)
      {
        break;
      }
    }
  }
  return true;
}
