#include "cellml/03_cellml_adapter.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/structured_regular_fixed.h"
#include "data_management/solution_vector_mapping.h"
#include "mesh/mesh_manager.h"

//#include <libcellml>    // libcellml not used here

template<int nStates, typename FunctionSpaceType>
CellmlAdapter<nStates,FunctionSpaceType>::
CellmlAdapter(DihuContext context) :
  CallbackHandler<nStates,FunctionSpaceType>(context)
{
  LOG(TRACE) << "CellmlAdapter constructor";
}

template<int nStates, typename FunctionSpaceType>
void CellmlAdapter<nStates,FunctionSpaceType>::
reset()
{
}
  
template<int nStates, typename FunctionSpaceType>
void CellmlAdapter<nStates,FunctionSpaceType>::
initialize()
{
  LOG(TRACE) << "CellmlAdapter<nStates,FunctionSpaceType>::initialize";

  CellmlAdapterBase<nStates,FunctionSpaceType>::initialize();
  
  // load rhs routine
  this->initializeRhsRoutine();

  this->initializeCallbackFunctions();
  
  int outputStateIndex = PythonUtility::getOptionInt(this->specificSettings_, "outputStateIndex", 0, PythonUtility::NonNegative);
  double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);

  // The solutionVectorMapping_ object stores the information which component of the solution will be further used
  // in methods that use the result of this method, e.g. in operator splitting.
  // The solutionVectorMapping object also scales the solution after transfer.
  this->solutionVectorMapping_.setOutputComponentNo(outputStateIndex);
  this->solutionVectorMapping_.setScalingFactor(prefactor);
}

template<int nStates, typename FunctionSpaceType>
void CellmlAdapter<nStates,FunctionSpaceType>::
initialize(double timeStepWidth)
{
}

template<int nStates, typename FunctionSpaceType>
void CellmlAdapter<nStates,FunctionSpaceType>::
setRankSubset(Partition::RankSubset rankSubset)
{
  // do nothing because we don't have stored data here (the data on which the computation is performed comes in evaluateTimesteppingRightHandSide from parameters) 
}

template<int nStates, typename FunctionSpaceType>
void CellmlAdapter<nStates,FunctionSpaceType>::
evaluateTimesteppingRightHandSideExplicit(Vec& input, Vec& output, int timeStepNo, double currentTime)
{
  //PetscUtility::getVectorEntries(input, states_);
  double *states, *rates;
  PetscErrorCode ierr;
  ierr = VecGetArray(input, &states); CHKERRV(ierr);   // get r/w pointer to contiguous array of the data, VecRestoreArray() needs to be called afterwards
  ierr = VecGetArray(output, &rates); CHKERRV(ierr);

  int nStatesInput, nRates;
  ierr = VecGetSize(input, &nStatesInput); CHKERRV(ierr);
  ierr = VecGetSize(output, &nRates); CHKERRV(ierr);
  VLOG(1) << "evaluateTimesteppingRightHandSideExplicit, input nStates: " << nStatesInput << ", output nRates: " << nRates;
  if (nStatesInput != nStates*this->nInstances_)
  {
    LOG(ERROR) << "nStatesInput=" << nStatesInput << ", nStates=" << nStates << ", nInstances=" << this->nInstances_;
  }
  assert (nStatesInput == nStates*this->nInstances_);
  assert (nRates == nStates*this->nInstances_);

  //LOG(DEBUG) << " evaluateTimesteppingRightHandSide: nInstances=" << this->nInstances_ << ", nStates=" << nStates;
  
  // get new values for parameters, call callback function of python config
  if (this->setParameters_ && timeStepNo % this->setParametersCallInterval_ == 0)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;
    
    VLOG(1) << "call setParameters";
    this->setParameters_((void *)this, this->nInstances_, timeStepNo, currentTime, this->parameters_);
  }

  //              this          STATES, RATES, WANTED,                KNOWN
  if(this->rhsRoutine_)
  {
    VLOG(1) << "call rhsRoutine_ with " << this->intermediates_.size() << " intermediates, " << this->parameters_.size() << " parameters";
    VLOG(2) << "intermediates: " << this->intermediates_ << ", parameters: " << this->parameters_;

    // call actual rhs routine from cellml code
    this->rhsRoutine_((void *)this, currentTime, states, rates, this->intermediates_.data(), this->parameters_.data());
  }

  // handle intermediates, call callback function of python config
  if (this->handleResult_ && timeStepNo % this->handleResultCallInterval_ == 0)
  {
    int nStatesInput;
    VecGetSize(input, &nStatesInput);

    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;
    
    VLOG(1) << "call handleResult with in total " << nStatesInput << " states, " << this->intermediates_.size() << " intermediates";    
    this->handleResult_((void *)this, this->nInstances_, timeStepNo, currentTime, states, this->intermediates_.data());
  }

  //PetscUtility::setVector(rates_, output);
  // give control of data back to Petsc
  ierr = VecRestoreArray(input, &states); CHKERRV(ierr);
  ierr = VecRestoreArray(output, &rates); CHKERRV(ierr);
}

/*
template<int nStates, typename FunctionSpaceType>
void CellmlAdapter<nStates,FunctionSpaceType>::
evaluateTimesteppingRightHandSideImplicit(Vec& input, Vec& output, int timeStepNo, double currentTime)
{
}
*/

//! return false because the object is independent of mesh type
template<int nStates, typename FunctionSpaceType>
bool CellmlAdapter<nStates,FunctionSpaceType>::
knowsMeshType()
{
  return CellmlAdapterBase<nStates,FunctionSpaceType>::knowsMeshType();
}

template<int nStates, typename FunctionSpaceType>
void CellmlAdapter<nStates,FunctionSpaceType>::
getComponentNames(std::vector<std::string> &stateNames)
{
  this->getStateNames(stateNames);
}

//! return the mesh
template<int nStates, typename FunctionSpaceType>
std::shared_ptr<Mesh::Mesh> CellmlAdapter<nStates,FunctionSpaceType>::
mesh()
{
  return CellmlAdapterBase<nStates,FunctionSpaceType>::mesh();
}

template<int nStates, typename FunctionSpaceType>
template<typename FunctionSpaceType2>
bool CellmlAdapter<nStates,FunctionSpaceType>::
setInitialValues(FieldVariable::FieldVariable<FunctionSpaceType2,nStates> &initialValues)
{
  return CellmlAdapterBase<nStates,FunctionSpaceType>::template setInitialValues<FunctionSpaceType2>(initialValues);
}
