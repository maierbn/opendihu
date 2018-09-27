#include "cellml/00_cellml_adapter_base.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/mesh_manager.h"

template<int nStates, typename FunctionSpaceType>
CellmlAdapterBase<nStates,FunctionSpaceType>::
CellmlAdapterBase(DihuContext context) :
  context_(context)
{
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "CellML");
  outputWriterManager_.initialize(specificSettings_);
  LOG(TRACE) << "CellmlAdapterBase constructor";
}

template<int nStates, typename FunctionSpaceType>
CellmlAdapterBase<nStates,FunctionSpaceType>::
~CellmlAdapterBase()
{
}

template<int nStates, typename FunctionSpaceType>
constexpr int CellmlAdapterBase<nStates,FunctionSpaceType>::
nComponents()
{
  return nStates;
}



template<int nStates, typename FunctionSpaceType>
void CellmlAdapterBase<nStates,FunctionSpaceType>::
initialize()
{
  LOG(TRACE) << "CellmlAdapterBase<nStates,FunctionSpaceType>::initialize";

  if (VLOG_IS_ON(1))
  {
    LOG(DEBUG) << "CellmlAdapterBase<nStates,FunctionSpaceType>::initialize querying meshManager for mesh";
    LOG(DEBUG) << "specificSettings_: ";
    PythonUtility::printDict(specificSettings_);
  }
  
  // create a mesh if there is not yet one assigned, function space FunctionSpace::Generic, downcasted to Mesh::Mesh
  functionSpace_ = context_.meshManager()->functionSpace<FunctionSpaceType>(specificSettings_);  // create initialized mesh
  LOG(DEBUG) << "Cellml mesh has " << functionSpace_->nNodesLocalWithoutGhosts() << " local nodes";

  //store number of instances
  nInstances_ = functionSpace_->nNodesLocalWithoutGhosts();

  stateNames_.resize(nStates);
  sourceFilename_ = PythonUtility::getOptionString(this->specificSettings_, "sourceFilename", "");
  this->scanSourceFile(this->sourceFilename_, statesInitialValues_);
  
  // add explicitely defined parameters that replace intermediates and constants
  if (!inputFileTypeOpenCMISS_)
  {
    PythonUtility::getOptionVector(this->specificSettings_, "parametersUsedAsIntermediate", parametersUsedAsIntermediate_);
    PythonUtility::getOptionVector(this->specificSettings_, "parametersUsedAsConstant", parametersUsedAsConstant_);
    nParameters_ += parametersUsedAsIntermediate_.size() + parametersUsedAsConstant_.size();
    
    //LOG(DEBUG) << "parametersUsedAsIntermediate_: " << parametersUsedAsIntermediate_ 
    //  << ", parametersUsedAsConstant_: " << parametersUsedAsConstant_;
  }
  
  LOG(DEBUG) << "Initialize CellML with nInstances = " << nInstances_ << ", nParameters_ = " << nParameters_ 
    << ", nStates = " << nStates << ", nIntermediates = " << nIntermediates_;
    
  // allocate data vectors
  intermediates_.resize(nIntermediates_*nInstances_);
  parameters_.resize(nParameters_*nInstances_);
  LOG(DEBUG) << "parameters.size: " << parameters_.size() << ", intermediates.size: " << intermediates_.size();
}


template<int nStates, typename FunctionSpaceType>
template<typename FunctionSpaceType2>
bool CellmlAdapterBase<nStates,FunctionSpaceType>::
setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nStates>> initialValues)
{
  LOG(TRACE) << "CellmlAdapterBase<nStates,FunctionSpaceType>::setInitialValues, sourceFilename_=" << this->sourceFilename_;
  if(PythonUtility::hasKey(this->specificSettings_, "statesInitialValues"))
  {
    LOG(DEBUG) << "set initial values from config";

    // statesInitialValues gives the initial state values for one instance of the problem. it is used for all instances.
    statesInitialValues_ = PythonUtility::getOptionArray<double,nStates>(this->specificSettings_, "statesInitialValues", 0);
  }
  else if(this->sourceFilename_ != "")
  {
    LOG(DEBUG) << "set initial values from source file";
    // parsing the source file was already done, the initial values are stored in the statesInitialValues_ vector
  }
  else
  {
    LOG(DEBUG) << "initialize to zero";
    statesInitialValues_.fill(0.0);
  }

  if(PythonUtility::hasKey(specificSettings_, "parametersInitialValues"))
  {
    LOG(DEBUG) << "load parametersInitialValues from config";

    std::vector<double> parametersInitial;
    PythonUtility::getOptionVector(specificSettings_, "parametersInitialValues", parametersInitial);

    VLOG(1) << "parametersInitialValues: " << parametersInitial << ", parameters_.size(): " << parameters_.size();

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
    VLOG(1) << "parameters_: " << parameters_;
  }
  else
  {
    LOG(DEBUG) << "Config does not contain key \"parametersInitialValues\"";
  }


  // Here we have the initial values for the states in the statesInitialValues_ vector, only for one instance.
  VLOG(1) << "statesInitialValues_: " << statesInitialValues_;

  const std::vector<std::array<double,nStates>> statesAllInstances(nInstances_, statesInitialValues_);

  VLOG(1) << "statesAllInstances: " << statesAllInstances << ", nInstances: " << nInstances_ << ", nStates per instances: " << statesInitialValues_.size();

  initialValues->setValuesWithoutGhosts(statesAllInstances);

  VLOG(1) << "initialValues: " << *initialValues;
  return true;

  LOG(DEBUG) << "do not set initial values";
  return false;
}

template<int nStates, typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> CellmlAdapterBase<nStates,FunctionSpaceType>::
functionSpace()
{
  return functionSpace_;
}

template<int nStates, typename FunctionSpaceType>
void CellmlAdapterBase<nStates,FunctionSpaceType>::
getNumbers(int& nInstances, int& nIntermediates, int& nParameters)
{
  nInstances = nInstances_;
  nIntermediates = nIntermediates_;
  nParameters = nParameters_;
}

template<int nStates, typename FunctionSpaceType>
void CellmlAdapterBase<nStates,FunctionSpaceType>::
getStateNames(std::vector<std::string> &stateNames)
{
  stateNames = this->stateNames_;
}

template<int nStates, typename FunctionSpaceType>
bool CellmlAdapterBase<nStates,FunctionSpaceType>::
knowsMeshType()
{
  return false;
}

