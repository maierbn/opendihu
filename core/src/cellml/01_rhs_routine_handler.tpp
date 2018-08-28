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
        doCompilation = false;
        file.close();
      }
    }
    
    if (doCompilation)
    {
     
      if (libraryFilename.find("/") != std::string::npos)
      {
        std::string path = libraryFilename.substr(0, libraryFilename.find("/"));
        int ret = system((std::string("mkdir -p ")+path).c_str());
        
        if (ret != 0)
        {
          LOG(ERROR) << "Could not create path \"" << path << "\".";
        }
      }
     
      std::stringstream compileCommand;
      // -ftree-vectorize -fopt-info-vec-missed -fopt-info-vec-optimized
      compileCommand << "gcc -fPIC -O3 -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared -lm -x c -o " << libraryFilename << " " << simdSourceFilename;
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
    rhsRoutine_ = (void (*)(double,double*,double*,double*,double*)) dlsym(handle, "OC_CellML_RHS_routine");
    rhsRoutineSimd_ = (void (*)(void *,double*,double*,double*,double*)) dlsym(handle, "OC_CellML_RHS_routine_simd");

    LOG(DEBUG) << "Library \"" <<libraryFilename<< "\" loaded. "
      << "rhsRoutine: " << (rhsRoutine_==NULL? "NULL" : "yes") << ", rhsRoutineSimd: " << (rhsRoutineSimd_==NULL? "NULL" : "yes");

    // fail if none of both could be loaded
    if (!rhsRoutine_ && !rhsRoutineSimd_)
    {
      LOG(FATAL) << "Could not load rhs routine from dynamic library \"" <<libraryFilename<< "\".";
    }

    // if only non-simd version could be loaded, create other from that
    if(!rhsRoutineSimd_)
    {
      LOG(DEBUG) << "Only non-simd rhs routine available";

      rhsRoutineSimd_ = [](void *context, double* OC_STATE, double* OC_RATE, double* OC_WANTED, double* OC_KNOWN)
      {
        LOG(DEBUG) << "call rhsRoutine ";

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
            states[i] = OC_STATE[i*nInstances + instanceNo];
          }
          for (int i = 0; i < nParameters; i++)
          {
            parameters[i] = OC_KNOWN[i*nInstances + instanceNo];
          }

          // call old rhs routine
          cellmlAdapter->rhsRoutine_(0.0, states.data(), rates.data(), intermediates.data(), parameters.data());

          // copy the results back in the big vectors

          // copy the values back to the interleaved vectors
          for (int i = 0; i < nStates; i++)
          {
            OC_RATE[i*nInstances + instanceNo] = rates[i];
          }

          for (int i = 0; i < nIntermediates; i++)
          {
            OC_WANTED[i*nInstances + instanceNo] = intermediates[i];
          }
        }
      };
    }
  }
  else
  {
    LOG(FATAL) << "Could not load dynamic library \"" <<libraryFilename<< "\". Reason: " <<dlerror();
  }
}

template<int nStates>
bool RhsRoutineHandler<nStates>::
createSimdSourceFile(std::string &simdSourceFilename)
{
  this->sourceFilename_ = PythonUtility::getOptionString(this->specificSettings_, "sourceFilename", "");

  // read in source from file
  std::ifstream sourceFile(this->sourceFilename_.c_str());
  if (!sourceFile.is_open())
  {
    LOG(ERROR) << "Could not open source file \"" <<this->sourceFilename_<< "\" for reading!";
    return false;
  }
  else
  {
    std::stringstream source;
    source << sourceFile.rdbuf();
    sourceFile.close();

    std::stringstream simdSource;

    // step through lines and create simd source
    while(!source.eof())
    {
      std::string line;
      getline(source, line);

      // line contains declaration of ALGEBRAIC variable
      if (line.find("double CONSTANTS") == 0 && line.find("ALGEBRAIC[") != std::string::npos)
      {
        size_t pos = line.find("ALGEBRAIC[");
        size_t posStart = line.find("[", pos)+1;
        size_t posEnd = line.find("]", posStart);
        int algebraicSize = atoi(line.substr(posStart).c_str());

        simdSource << line.substr(0, posStart) << algebraicSize * this->nInstances_ << line.substr(posEnd) << std::endl;
      }
      // line contains function head
      else if (line.find("void OC_CellML_RHS_routine") != std::string::npos)
      {
        if (line.find("void OC_CellML_RHS_routine_simd") != std::string::npos)
        {
          LOG(WARNING) << "The given source file \"" <<this->sourceFilename_<< "\" already contains a simd version of the rhs routine. "
            << "Use the option \"simdSourceFilename\" instead.";
          return true;
        }

        simdSource << "/* Routine is designed for " << this->nInstances_ << " instances of the CellML problem. */" << std::endl
          << "void OC_CellML_RHS_routine_simd("
          << "void *context, double* OC_STATE, double* OC_RATE, double* OC_WANTED, double* OC_KNOWN)" << std::endl;
      }
      // line contains assignment
      else if(line.find("OC_WANTED") == 0 || line.find("OC_RATE") == 0 || line.find("ALGEBRAIC") == 0)
      {
        // parse line
        struct entry_t
        {
          std::string code;
          int arrayIndex;
          enum {variableName, other} type;
        };
        std::list<entry_t> entries;

        VLOG(2) << "line: [" <<line<< "]";

        size_t currentPos = 0;
        while(currentPos <= line.length())
        {
          VLOG(2);
          VLOG(2) << "currentPos: " << currentPos<< " code from there: \"" <<line.substr(currentPos, 20) << "\"";
          VLOG(2) << "variables: "
            <<line.find("OC_STATE", currentPos) << ", "
            <<line.find("OC_RATE", currentPos) << ", "
            <<line.find("ALGEBRAIC", currentPos) << ", "
            <<line.find("OC_WANTED", currentPos) << ", "
            <<line.find("OC_KNOWN", currentPos);

          size_t posVariable = std::min({
            line.find("OC_STATE", currentPos),
            line.find("OC_RATE", currentPos),
            line.find("ALGEBRAIC", currentPos),
            line.find("OC_WANTED", currentPos),
            line.find("OC_KNOWN", currentPos)
          });

          VLOG(2) << "posVariable: " <<posVariable;

          entry_t entry;
          if (posVariable == currentPos)
          {
            size_t posBracket = line.find("[", currentPos);
            entry.code = line.substr(currentPos, posBracket-currentPos);
            entry.arrayIndex = atoi(line.substr(posBracket+1).c_str());

            VLOG(2) << "extract variable name \"" <<entry.code<< "\", index " <<entry.arrayIndex;

            entry.type = entry_t::variableName;

            if (entry.code == "OC_STATE" && entry.arrayIndex >= nStates)
            {
              LOG(FATAL) << "CellML code in source file \"" << this->sourceFilename_ << "\" "
                << "computes more states than given in the user code as template parameter "
                << "(template parameter: " << nStates << ", CellML code: at least " << entry.arrayIndex+1 << ").";
            }

            // advance current position
            currentPos = line.find("]", currentPos)+1;
            VLOG(2) << "(1)advance to " << currentPos;
          }
          else
          {
            entry.code = line.substr(currentPos, posVariable-currentPos);
            VLOG(2) << "extract code \"" <<entry.code<< "\".";
            entry.type = entry_t::other;

            // advance current position
            currentPos = posVariable;
            VLOG(2) << "(2)advance to " << currentPos;
          }
          entries.push_back(entry);
        }

        simdSource << "for(int i = 0; i < " << this->nInstances_ << "; i++)" << std::endl
          << "{" << std::endl << "\t";

        VLOG(2) << "parsed " <<entries.size() << " entries";

        // write out parsed code with adjusted indices
        for (auto entry : entries)
        {
          VLOG(2) << " entry type=" <<entry.type<< ", code: \"" <<entry.code<< "\", index: " <<entry.arrayIndex;
          switch(entry.type)
          {
            case entry_t::variableName:
              simdSource << entry.code << "[" << entry.arrayIndex * this->nInstances_ << "+i]";
              break;
            case entry_t::other:
              simdSource << entry.code;
              break;
          }
        }
        simdSource << std::endl << "}" << std::endl;
      }
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
      LOG(ERROR) << "Could not write to file \"" <<simdSourceFilename<< "\".";
      return false;
    }
    else
    {
      simdSourceFile << simdSource.str();
      simdSourceFile.close();
    }
  }
  
  return true;
}