#include "cellml/cellml_adapter.h"

#include <list>

#include <Python.h>

#include "control/python_utility.h"
#include "control/petsc_utility.h"
#include "mesh/regular_fixed.h"
#include "data_management/solution_vector_mapping.h"
#include "mesh/mesh_manager.h"

//#include <libcellml>    // libcellml not used here

#include <unistd.h>
#include <dlfcn.h>

CellmlAdapter::CellmlAdapter(const DihuContext& context) :
  DiscretizableInTime(SolutionVectorMapping(true)),
  context_(context), setParameters_(NULL), handleResult_(NULL), 
  pythonSetParametersFunction_(NULL), pythonHandleResultFunction_(NULL)
{
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "CellML");
  outputWriterManager_.initialize(specificSettings_);
}

CellmlAdapter::~CellmlAdapter()
{
  Py_CLEAR(pythonSetParametersFunction_);
  Py_CLEAR(pythonHandleResultFunction_);
}


void CellmlAdapter::registerHandleResult(void (*handleResult) (void *context, int nInstances, int timeStepNo, double currentTime, 
                                                                 double *states, double *intermediates))
{
  handleResult_ = handleResult;
}

void CellmlAdapter::registerSetParameters(void (*setParameters) (void *context, int nInstances, int timeStepNo, double currentTime, 
                                                                 std::vector<double> &parmeters))
{
  setParameters_ = setParameters;
}

bool CellmlAdapter::createSimdSourceFile(std::string &simdSourceFilename)
{
  std::string sourceFilename = PythonUtility::getOptionString(specificSettings_, "sourceFilename", "");
  
  // read in source from file
  std::ifstream sourceFile(sourceFilename.c_str());
  if (!sourceFile.is_open())
  {
    LOG(ERROR) << "Could not open source file \""<<sourceFilename<<"\" for reading!";
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
        
        simdSource << line.substr(0, posStart) << algebraicSize * nInstances_ << line.substr(posEnd) << std::endl;
      }
      // line contains function head
      else if (line.find("void OC_CellML_RHS_routine") != std::string::npos)
      {
        if (line.find("void OC_CellML_RHS_routine_simd") != std::string::npos)
        {
          LOG(WARNING) << "The given source file \""<<sourceFilename<<"\" already contains a simd version of the rhs routine. "
            << "Use the option \"simdSourceFilename\" instead.";
          return true;
        }
        
        simdSource << "/* Routine is designed for " << nInstances_ << " instances of the CellML problem. */" << std::endl
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
        
        VLOG(2) << "line: ["<<line<<"]";
        
        size_t currentPos = 0;
        while(currentPos <= line.length())
        {
          VLOG(2);
          VLOG(2) << "currentPos: "<<currentPos<<" code from there: \""<<line.substr(currentPos, 20)<<"\"";
          VLOG(2) << "variables: "
            <<line.find("OC_STATE", currentPos)<<", "
            <<line.find("OC_RATE", currentPos)<<", "
            <<line.find("ALGEBRAIC", currentPos)<<", "
            <<line.find("OC_WANTED", currentPos)<<", "
            <<line.find("OC_KNOWN", currentPos);
          
          size_t posVariable = std::min({
            line.find("OC_STATE", currentPos), 
            line.find("OC_RATE", currentPos),
            line.find("ALGEBRAIC", currentPos), 
            line.find("OC_WANTED", currentPos), 
            line.find("OC_KNOWN", currentPos)
          });
          
          VLOG(2) << "posVariable: "<<posVariable;
          
          entry_t entry;
          if (posVariable == currentPos)
          {
            size_t posBracket = line.find("[", currentPos);
            entry.code = line.substr(currentPos, posBracket-currentPos);
            entry.arrayIndex = atoi(line.substr(posBracket+1).c_str());
            
            VLOG(2) << "extract variable name \""<<entry.code<<"\", index "<<entry.arrayIndex;
            
            entry.type = entry_t::variableName;
            
            if (entry.code == "OC_STATE" && entry.arrayIndex >= nStates_)
            {
              LOG(FATAL) << "CellML code in source file \"" << sourceFilename << "\" " 
                << "computes more states than specified in config " 
                << "(config: " << nStates_ << ", code: at least " << entry.arrayIndex+1 << ").";
            }
            
            // advance current position
            currentPos = line.find("]", currentPos)+1;
            VLOG(2) << "(1)advance to "<<currentPos;
          }
          else
          {
            entry.code = line.substr(currentPos, posVariable-currentPos);
            VLOG(2) << "extract code \""<<entry.code<<"\".";
            entry.type = entry_t::other;
            
            // advance current position
            currentPos = posVariable;
            VLOG(2) << "(2)advance to "<<currentPos;
          }
          entries.push_back(entry);
        }
        
        simdSource << "for(int i = 0; i < " << nInstances_ << "; i++)" << std::endl
          << "{" << std::endl << "\t";
          
        VLOG(2) << "parsed "<<entries.size()<<" entries";
          
        // write out parsed code with adjusted indices
        for (auto entry : entries)
        {
          VLOG(2) << " entry type=" <<entry.type<<", code: \""<<entry.code<<"\", index: "<<entry.arrayIndex;
          switch(entry.type)
          {
            case entry_t::variableName:
              simdSource << entry.code << "[" << entry.arrayIndex * nInstances_ << "+i]";
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
    simdSourceFilename = "simd_source.c";
    if (PythonUtility::containsKey(specificSettings_, "simdSourceFilename"))
    {
      simdSourceFilename = PythonUtility::getOptionString(specificSettings_, "simdSourceFilename", "");
    }
    
    std::ofstream simdSourceFile(simdSourceFilename.c_str());
    if (!simdSourceFile.is_open())
    {
      LOG(ERROR) << "Could not write to file \""<<simdSourceFilename<<"\".";
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

bool CellmlAdapter::scanInitialValues(std::string sourceFilename, std::vector<double> &statesInitialValues)
{
  LOG(TRACE) << "scanInitialValues";
  
  // read in source from file
  std::ifstream sourceFile(sourceFilename.c_str());
  if (!sourceFile.is_open())
  {
    LOG(WARNING) << "Could not open source file \""<<sourceFilename<<"\" for reading initial values.";
    return false;
  }
  else 
  {
    std::stringstream source;
    source << sourceFile.rdbuf();
    sourceFile.close();
    
    statesInitialValues.resize(nStates_);
    
    // step through lines and create simd source
    while(!source.eof())
    {
      std::string line;
      getline(source, line);
      
      // line contains initial value for a state or a known value, for example: "DUMMY_ASSIGNMENT /*OC_STATE[0]*/ = -79.974;"
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
        
        if (variableName == "OC_STATE" && index >= 0 && index < (unsigned int)nStates_)
        {
          statesInitialValues[index] = value;
        }
        else if (variableName == "OC_KNOWN" && index >= 0)
        {
          if (parameters_.size() < index+1)
            parameters_.resize(index+1);
          parameters_[index] = value;
        }
      }
    }
  }
  return true;
}

void CellmlAdapter::initializeRhsRoutine()
{
  // 1) if simdSourceFilename is given, use that source to compile the library
  // 2) if not 1) but sourceFilename is given, create simdSourceFilename from that and compile library
  // 3) if not 2) but libraryFilename is given, load that library, if it contains simdRhs, use that, if it contains non-simd rhs use that
  
  // try to load or create simd source
  std::string simdSourceFilename;
  if (PythonUtility::containsKey(specificSettings_, "sourceFilename"))
  {
    if (!createSimdSourceFile(simdSourceFilename))
      simdSourceFilename = "";
  }
  
  if (simdSourceFilename == "")
  {
    if(PythonUtility::containsKey(specificSettings_, "simdSourceFilename"))
    {
      simdSourceFilename = PythonUtility::getOptionString(specificSettings_, "simdSourceFilename", "");
    }
  }
  
  std::string libraryFilename;
  // if simdSourceFilename is set, compile to create dynamic library
  if (simdSourceFilename != "")
  {
    sourceFilename_ = simdSourceFilename;
    
    // compile source file to a library 
    libraryFilename = "lib.so";
    if (PythonUtility::containsKey(specificSettings_, "libraryFilename"))
    {
      libraryFilename = PythonUtility::getOptionString(specificSettings_, "libraryFilename", "lib.so");
    }
    
    std::stringstream compileCommand;
    // -ftree-vectorize -fopt-info-vec-missed -fopt-info-vec-optimized
    compileCommand << "gcc -fPIC -O3 -ftree-vectorize -fopt-info-vec-all=vectorizer_all.log -shared -lm -x c -o " << libraryFilename << " " << simdSourceFilename;
    int ret = system(compileCommand.str().c_str());
    if (ret != 0)
    {
      LOG(ERROR) << "Compilation failed. Command: \""<<compileCommand.str()<<"\".";
      libraryFilename = "";
    }
    else
    {
      LOG(DEBUG) << "Compilation successful. Command: \""<<compileCommand.str()<<"\".";
    }

    // repeat compilation with different GCC vectorizer outputs
#if 0    
    compileCommand.str("");
    compileCommand << "gcc -fPIC -O3 -ftree-vectorize -fopt-info-vec-missed=vectorizer_missed.log -shared -lm -x c -o " << libraryFilename << " " << simdSourceFilename;
    ret = system(compileCommand.str().c_str());
    if (ret != 0)
    {
      LOG(DEBUG) << "Compilation failed. Command: \""<<compileCommand.str()<<"\".";
    }
    else
    {
      LOG(DEBUG) << "Compilation successful. Command: \""<<compileCommand.str()<<"\".";
    }
    compileCommand.str("");
    compileCommand << "gcc -fPIC -O3 -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared -lm -x c -o " << libraryFilename << " " << simdSourceFilename;
    ret = system(compileCommand.str().c_str());
    if (ret != 0)
    {
      LOG(DEBUG) << "Compilation failed. Command: \""<<compileCommand.str()<<"\".";
    }
    else
    {
      LOG(DEBUG) << "Compilation successful. Command: \""<<compileCommand.str()<<"\".";
    }
#endif
  }
  
  // if library is still not available, look if it was provided
  if (libraryFilename == "")
  {
    libraryFilename = PythonUtility::getOptionString(specificSettings_, "libraryFilename", "lib.so");
  }
  
  if (libraryFilename == "")
  {
    LOG(FATAL) << "Could not create or locate dynamic library for cellml right hand side routine.";
  }
  
  // load dynamic library
  void* handle = dlopen(libraryFilename.c_str(), RTLD_LOCAL | RTLD_LAZY);
  if (handle)
  {
    // load rhs method
    rhsRoutine_ = (void (*)(double,double*,double*,double*,double*)) dlsym(handle, "OC_CellML_RHS_routine");
    rhsRoutineSimd_ = (void (*)(void *,double*,double*,double*,double*)) dlsym(handle, "OC_CellML_RHS_routine_simd");
   
    LOG(DEBUG) << "Library \""<<libraryFilename<<"\" loaded. " 
      << "rhsRoutine: " << (rhsRoutine_==NULL? "NULL" : "yes") << ", rhsRoutineSimd: " << (rhsRoutineSimd_==NULL? "NULL" : "yes");
    
    // fail if none of both could be loaded
    if (!rhsRoutine_ && !rhsRoutineSimd_)
    {
      LOG(FATAL) << "Could not load rhs routine from dynamic library \""<<libraryFilename<<"\".";
    }
    
    // if only non-simd version could be loaded, create other from that
    if(!rhsRoutineSimd_)
    {
      LOG(DEBUG) << "Only non-simd rhs routine available";
      
      rhsRoutineSimd_ = [](void *context, double* OC_STATE, double* OC_RATE, double* OC_WANTED, double* OC_KNOWN)
      {
        LOG(DEBUG) << "call rhsRoutine ";
        
        CellmlAdapter *cellmlAdapter = (CellmlAdapter *)context;
        int nStates, nInstances, nIntermediates, nParameters;
        cellmlAdapter->getNumbers(nStates, nInstances, nIntermediates, nParameters);
        
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
    LOG(FATAL) << "Could not load dynamic library \""<<libraryFilename<<"\".";
  }
    
}

void CellmlAdapter::initialize()
{
  LOG(TRACE) << "CellmlAdapter::initialize";
  
  // parse number of variables
  nStates_ = PythonUtility::getOptionInt(specificSettings_, "numberStates", 0, PythonUtility::NonNegative);
  nIntermediates_ = PythonUtility::getOptionInt(specificSettings_, "numberIntermediates", 0, PythonUtility::NonNegative);
  nParameters_ = PythonUtility::getOptionInt(specificSettings_, "numberParameters", 0, PythonUtility::NonNegative);
  
  LOG(DEBUG) << "CellmlAdapter::initialize querying meshManager for mesh";
  LOG(DEBUG) << "specificSettings_: ";
  PythonUtility::printDict(specificSettings_);
  
  // create a mesh if there is not yet one assigned
  mesh_ = context_.meshManager()->mesh<Mesh::Mesh>(specificSettings_);
  LOG(DEBUG) << "Cellml mesh has " << mesh_->nNodes() << " nodes";
  
  //store number of instances
  nInstances_ = mesh_->nNodes();
  
  LOG(DEBUG) << "Initialize CellML with nStates="<<nStates_
    <<", nIntermediates="<<nIntermediates_<<", nParameters="<<nParameters_<<", nInstances="<<nInstances_;
  
  // load rhs routine
  initializeRhsRoutine();
  
  LOG(DEBUG) << "rhs initialized";
  
  // allocate data vectors
  intermediates_.resize(nIntermediates_*nInstances_);
  parameters_.resize(nParameters_*nInstances_);
  LOG(DEBUG) << "size of parameters: "<<parameters_.size();
  
  int outputStateIndex = PythonUtility::getOptionInt(specificSettings_, "outputStateIndex", 0, PythonUtility::NonNegative);
  double prefactor = PythonUtility::getOptionDouble(specificSettings_, "prefactor", 1.0);
  
  // The solutionVectorMapping_ object stores the information which range of values of the solution will be further used 
  // in methods that use the result of this method, e.g. in operator splittings. 
  // These are all values of a single STATE with number outputStateIndex from settings.
  // The data layout is for e.g. 3 instances like this: STATE[0] STATE[0] STATE[0] STATE[1] STATE[1] STATE[1] STATE[2]...
  solutionVectorMapping_.setOutputRange(nInstances_*outputStateIndex, nInstances_*(outputStateIndex+1));
  solutionVectorMapping_.setScalingFactor(prefactor);
  

  if (PythonUtility::containsKey(specificSettings_, "setParametersFunction"))
  {
    pythonSetParametersFunction_ = PythonUtility::getOptionFunction(specificSettings_, "setParametersFunction");
    setParametersCallInterval_ = PythonUtility::getOptionInt(specificSettings_, "setParametersCallInterval", 1, PythonUtility::Positive);
    setParameters_ = [](void *context, int nInstances, int timeStepNo, double currentTime, std::vector<double> &parameters)
    {
      CellmlAdapter *cellmlAdapter = (CellmlAdapter *)context;
      cellmlAdapter->callPythonSetParametersFunction(nInstances, timeStepNo, currentTime, parameters);
    };
    LOG(DEBUG) << "registered setParameters function";
  }

  if (PythonUtility::containsKey(specificSettings_, "handleResultFunction"))
  {
    pythonHandleResultFunction_ = PythonUtility::getOptionFunction(specificSettings_, "handleResultFunction");
    handleResultCallInterval_ = PythonUtility::getOptionInt(specificSettings_, "handleResultCallInterval", 1, PythonUtility::Positive);
    handleResult_ = [](void *context, int nInstances, int timeStepNo, double currentTime, double *states, double *intermediates)
    {
      CellmlAdapter *cellmlAdapter = (CellmlAdapter *)context;
      cellmlAdapter->callPythonHandleResultFunction(nInstances, timeStepNo, currentTime, states, intermediates);
    };
    LOG(DEBUG) << "registered handleResult function";
  }
}

void CellmlAdapter::callPythonSetParametersFunction(int nInstances, int timeStepNo, double currentTime, std::vector< double >& parameters)
{
  if (pythonSetParametersFunction_ == NULL)
    return;
  
  // compose callback function
  PyObject *parametersList = PythonUtility::convertToPythonList(parameters);
  PyObject *arglist = Py_BuildValue("(i,i,d,O)", nInstances, timeStepNo, currentTime, parametersList);
  PyObject *returnValue = PyObject_CallObject(pythonSetParametersFunction_, arglist);
  
  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();
  
  // copy new values in parametersList to parameters_ vector
  for (unsigned int i=0; i<parameters.size(); i++)
  {
    PyObject *item = PyList_GetItem(parametersList, (Py_ssize_t)i);
    parameters[i] = PythonUtility::convertFromPython<double>(item);
  }
  
  // decrement reference counters for python objects
  Py_CLEAR(parametersList);
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist); 
}

void CellmlAdapter::callPythonHandleResultFunction(int nInstances, int timeStepNo, double currentTime, 
                                                    double *states, double *intermediates)
{
  if (pythonHandleResultFunction_ == NULL)
    return;
  
  // compose callback function
  LOG(DEBUG) << "callPythonHandleResultFunction: nInstances: " << nInstances_<<", nStates: " << nStates_ << ", nIntermediates: " << nIntermediates_;
  PyObject *statesList = PythonUtility::convertToPythonList(nStates_*nInstances_, states);
  PyObject *intermediatesList = PythonUtility::convertToPythonList(nIntermediates_*nInstances_, intermediates);
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O)", nInstances, timeStepNo, currentTime, statesList, intermediatesList);
  PyObject *returnValue = PyObject_CallObject(pythonHandleResultFunction_, arglist);
  
  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();
  
  // decrement reference counters for python objects
  Py_CLEAR(statesList);
  Py_CLEAR(intermediatesList);
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

bool CellmlAdapter::setInitialValues(Vec& initialValues)
{
  LOG(TRACE) << "CellmlAdapter::setInitialValues, sourceFilename_="<<sourceFilename_;
  std::vector<double> states;
  if(PythonUtility::containsKey(specificSettings_, "statesInitialValues"))
  {
    LOG(DEBUG) << "set initial values from config";

    PythonUtility::getOptionVector(specificSettings_, "statesInitialValues", nStates_, states);
  }
  else if(sourceFilename_ != "")
  {
    LOG(DEBUG) << "set initial values from source file";
    
    scanInitialValues(sourceFilename_, states);
  }
  else 
  {
    LOG(DEBUG) << "initialize to zero";
    states.resize(nStates_*nInstances_, 0);
  }
  
  if(PythonUtility::containsKey(specificSettings_, "parametersInitialValues"))
  {
    LOG(DEBUG) << "load parametersInitialValues also from config";
    
    std::vector<double> parametersInitial;
    PythonUtility::getOptionVector(specificSettings_, "parametersInitialValues", nStates_, parametersInitial);
    
    if (parametersInitial.size() == parameters_.size())
    {
      std::copy(parametersInitial.begin(), parametersInitial.end(), parameters_.begin());
      LOG(DEBUG) << "parameters size is matching for all instances";
    }
    else
    {
      LOG(DEBUG) << "copy parameters which were given only for one instance to all instances";
      for(int instanceNo=0; instanceNo<nInstances_; instanceNo++)
      {
        for(int j=0; j<nParameters_; j++)
        {
          parameters_[j*nInstances_ + instanceNo] = parametersInitial[j];
        }
      }
    }
  }
  else
  {
    LOG(DEBUG) << "Config does not contain key \"parametersInitialValues\"";
  }
  
  if (!states.empty())
  {
    std::vector<double> statesAllInstances(nStates_*nInstances_);
    for(int j=0; j<nStates_; j++)
    {
      for(int instanceNo=0; instanceNo<nInstances_; instanceNo++)
      {
        statesAllInstances[j*nInstances_ + instanceNo] = states[j];
      }
    }
    
    PetscUtility::setVector(statesAllInstances, initialValues);
      
    LOG(DEBUG) << "initial values were set as follows: ";
    for(auto value : statesAllInstances)
      LOG(DEBUG) << "  " << value;
    return true;
  }
  
  LOG(DEBUG) << "do not set initial values";
  
  return false;
}

std::shared_ptr<Mesh::Mesh> CellmlAdapter::mesh()
{
  return mesh_;
}

int CellmlAdapter::numberDegreesOfFreedomPerNode()
{
  // this is the number of entries per mesh node that the input and output vectors of evaluateTimesteppingRightHandSide will have
  return nStates_;
}

void CellmlAdapter::getNumbers(int& nStates, int& nInstances, int& nIntermediates, int& nParameters)
{
  nStates = nStates_;
  nInstances = nInstances_;
  nIntermediates = nIntermediates_;
  nParameters = nParameters_;
}

void CellmlAdapter::evaluateTimesteppingRightHandSide(Vec& input, Vec& output, int timeStepNo, double currentTime)
{
  //PetscUtility::getVectorEntries(input, states_);
  double *states, *rates;
  VecGetArray(input, &states);    // get r/w pointer to contiguous array of the data, VecRestoreArray() needs to be called afterwards
  VecGetArray(output, &rates);
  
  // get new values for parameters
  if (setParameters_ && timeStepNo % setParametersCallInterval_ == 0)
    setParameters_((void *)this, nInstances_, timeStepNo, currentTime, parameters_);
    
  //              this          STATES, RATES, WANTED,                KNOWN
  if(rhsRoutineSimd_)
    rhsRoutineSimd_((void *)this, states, rates, intermediates_.data(), parameters_.data());
  
  // handle intermediates
  if (handleResult_ && timeStepNo % handleResultCallInterval_ == 0)
  {
    int nStates;
    VecGetSize(input, &nStates);
    LOG(DEBUG) << "call handleResult with in total " << nStates << " states, " << intermediates_.size() << " intermediates";
    handleResult_((void *)this, nInstances_, timeStepNo, currentTime, states, intermediates_.data());
  }
  
  //PetscUtility::setVector(rates_, output);
  // give control of data back to Petsc
  VecRestoreArray(input, &states);
  VecRestoreArray(output, &rates);
}

bool CellmlAdapter::knowsMeshType()
{
  return false;
}
  