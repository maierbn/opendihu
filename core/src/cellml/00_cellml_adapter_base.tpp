#include "cellml/00_cellml_adapter_base.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/mesh_manager.h"

template<int nStates>
CellmlAdapterBase<nStates>::
CellmlAdapterBase(DihuContext context) :
  context_(context)
{
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "CellML");
  outputWriterManager_.initialize(specificSettings_);
  LOG(TRACE) << "CellmlAdapterBase constructor";
}

template<int nStates>
CellmlAdapterBase<nStates>::
~CellmlAdapterBase()
{
}

template<int nStates>
constexpr int CellmlAdapterBase<nStates>::
nComponents()
{
  return nStates;
}


template<int nStates>
bool CellmlAdapterBase<nStates>::
scanInitialValues(std::string sourceFilename, std::vector<double> &statesInitialValues)
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

    statesInitialValues.resize(nStates);

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

        if (variableName == "OC_STATE" && index >= 0 && index < (unsigned int)nStates)
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

template<int nStates>
void CellmlAdapterBase<nStates>::
initialize()
{
  LOG(TRACE) << "CellmlAdapterBase<nStates>::initialize";

  // parse number of variables
  nIntermediates_ = PythonUtility::getOptionInt(specificSettings_, "numberIntermediates", 0, PythonUtility::NonNegative);
  nParameters_ = PythonUtility::getOptionInt(specificSettings_, "numberParameters", 0, PythonUtility::NonNegative);
 
  if (VLOG_IS_ON(1))
  {
    LOG(DEBUG) << "CellmlAdapterBase<nStates>::initialize querying meshManager for mesh";
    LOG(DEBUG) << "specificSettings_: ";
    PythonUtility::printDict(specificSettings_);
  }
  
  // create a mesh if there is not yet one assigned
  mesh_ = context_.meshManager()->mesh<>(specificSettings_);
  mesh_->initialize();
  LOG(DEBUG) << "Cellml mesh has " << mesh_->nLocalNodes() << " local nodes";

  //store number of instances
  nInstances_ = mesh_->nLocalNodes();

  LOG(DEBUG) << "Initialize CellML with nStates="<<nStates
    <<", nIntermediates="<<nIntermediates_<<", nParameters="<<nParameters_<<", nInstances="<<nInstances_;

  // allocate data vectors
  intermediates_.resize(nIntermediates_*nInstances_);
  parameters_.resize(nParameters_*nInstances_);
  LOG(DEBUG) << "size of parameters: "<<parameters_.size();
}

template<int nStates>
bool CellmlAdapterBase<nStates>::
setInitialValues(Vec& initialValues)
{
  LOG(TRACE) << "CellmlAdapterBase<nStates>::setInitialValues, sourceFilename_="<<this->sourceFilename_;
  std::vector<double> states;
  if(PythonUtility::hasKey(this->specificSettings_, "statesInitialValues"))
  {
    LOG(DEBUG) << "set initial values from config";

    PythonUtility::getOptionVector(this->specificSettings_, "statesInitialValues", nStates, states);
  }
  else if(this->sourceFilename_ != "")
  {
    LOG(DEBUG) << "set initial values from source file";

    scanInitialValues(sourceFilename_, states);
  }
  else
  {
    LOG(DEBUG) << "initialize to zero";
    states.resize(nStates*nInstances_, 0);
  }

  if(PythonUtility::hasKey(specificSettings_, "parametersInitialValues"))
  {
    LOG(DEBUG) << "load parametersInitialValues also from config";

    std::vector<double> parametersInitial;
    PythonUtility::getOptionVector(specificSettings_, "parametersInitialValues", nStates, parametersInitial);

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
    std::vector<double> statesAllInstances(nStates*nInstances_);
    for(int j=0; j<nStates; j++)
    {
      for(int instanceNo=0; instanceNo<nInstances_; instanceNo++)
      {
        statesAllInstances[j*nInstances_ + instanceNo] = states[j];
      }
    }

    PetscUtility::setVector(statesAllInstances, initialValues);

    if (VLOG_IS_ON(2))
    {
      VLOG(2) << "initial values were set as follows: ";
      for(auto value : statesAllInstances)
        VLOG(2) << "  " << value;
    }
    return true;
  }

  LOG(DEBUG) << "do not set initial values";

  return false;
}

template<int nStates>
std::shared_ptr<Mesh::Mesh> CellmlAdapterBase<nStates>::
mesh()
{
  return mesh_;
}

template<int nStates>
void CellmlAdapterBase<nStates>::
getNumbers(int& nInstances, int& nIntermediates, int& nParameters)
{
  nInstances = nInstances_;
  nIntermediates = nIntermediates_;
  nParameters = nParameters_;
}

template<int nStates>
bool CellmlAdapterBase<nStates>::
knowsMeshType()
{
  return false;
}

