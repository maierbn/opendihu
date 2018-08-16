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

template<int nStates>
CellmlAdapter<nStates>::
CellmlAdapter(DihuContext context) :
  CallbackHandler<nStates>(context),
  DiscretizableInTime(SolutionVectorMapping(true))
{
  LOG(TRACE) << "CellmlAdapter constructor";
}

template<int nStates>
void CellmlAdapter<nStates>::
initialize()
{
  LOG(TRACE) << "CellmlAdapter<nStates>::initialize";

  CellmlAdapterBase<nStates>::initialize();
  
  // load rhs routine
  this->initializeRhsRoutine();

  this->initializeCallbackFunctions();
  
  int outputStateIndex = PythonUtility::getOptionInt(this->specificSettings_, "outputStateIndex", 0, PythonUtility::NonNegative);
  double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);

  // The solutionVectorMapping_ object stores the information which range of values of the solution will be further used
  // in methods that use the result of this method, e.g. in operator splittings.
  // These are all values of a single STATE with number outputStateIndex from settings.
  // The data layout is for e.g. 3 instances like this: STATE[0] STATE[0] STATE[0] STATE[1] STATE[1] STATE[1] STATE[2]...
  solutionVectorMapping_.setOutputRange(this->nInstances_*outputStateIndex, this->nInstances_*(outputStateIndex+1));
  solutionVectorMapping_.setScalingFactor(prefactor);
}

template<int nStates>
void CellmlAdapter<nStates>::
initialize(double timeStepWidth)
{
}


template<int nStates>
void CellmlAdapter<nStates>::
evaluateTimesteppingRightHandSideExplicit(Vec& input, Vec& output, int timeStepNo, double currentTime)
{
  //PetscUtility::getVectorEntries(input, states_);
  double *states, *rates;
  VecGetArray(input, &states);    // get r/w pointer to contiguous array of the data, VecRestoreArray() needs to be called afterwards
  VecGetArray(output, &rates);

  //LOG(DEBUG) << " evaluateTimesteppingRightHandSide: nInstances=" << this->nInstances_ << ", nStates=" << nStates;
  
  // get new values for parameters, call callback function of python config
  if (this->setParameters_ && timeStepNo % this->setParametersCallInterval_ == 0)
  {
    LOG(DEBUG) << "call setParameters";
    
    // start critical section for python interpreter (only one thread)
    //Py_BEGIN_ALLOW_THREADS
    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();
    
    this->setParameters_((void *)this, this->nInstances_, timeStepNo, currentTime, this->parameters_);
    
    
    /* Release the thread. No Python API allowed beyond this point. */
    PyGILState_Release(gstate);
    //Py_END_ALLOW_THREADS
  }

  //              this          STATES, RATES, WANTED,                KNOWN
  if(this->rhsRoutineSimd_)
  { 
    this->rhsRoutineSimd_((void *)this, states, rates, this->intermediates_.data(), this->parameters_.data());
  }

  // handle intermediates, call callback function of python config
  if (this->handleResult_ && timeStepNo % this->handleResultCallInterval_ == 0)
  {
    int nStatesInput;
    VecGetSize(input, &nStatesInput);

    LOG(DEBUG) << "call handleResult with in total " << nStatesInput << " states, " << this->intermediates_.size() << " intermediates";
    
    // start critical section for python interpreter (only one thread)
    //Py_BEGIN_ALLOW_THREADS
    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();
    
    this->handleResult_((void *)this, this->nInstances_, timeStepNo, currentTime, states, this->intermediates_.data());
    
    
    /* Release the thread. No Python API allowed beyond this point. */
    PyGILState_Release(gstate);
    //Py_END_ALLOW_THREADS
  }

  //PetscUtility::setVector(rates_, output);
  // give control of data back to Petsc
  VecRestoreArray(input, &states);
  VecRestoreArray(output, &rates);
}

template<int nStates>
void CellmlAdapter<nStates>::
evaluateTimesteppingRightHandSideImplicit(Vec& input, Vec& output, int timeStepNo, double currentTime)
{
}

//! return false because the object is independent of mesh type
template<int nStates>
bool CellmlAdapter<nStates>::
knowsMeshType()
{
  return CellmlAdapterBase<nStates>::knowsMeshType();
}

//! return the mesh

template<int nStates>
std::shared_ptr<Mesh::Mesh> CellmlAdapter<nStates>::
mesh()
{
  return CellmlAdapterBase<nStates>::mesh();
}

template<int nStates>
bool CellmlAdapter<nStates>::
setInitialValues(Vec &initialValues)
{
  return CellmlAdapterBase<nStates>::setInitialValues(initialValues);
}