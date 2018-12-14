#include "time_stepping_scheme/time_stepping_scheme_ode.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{
  template<typename FunctionSpaceType, int nComponents>
  TimeSteppingSchemeOdeBase<FunctionSpaceType,nComponents>::
  TimeSteppingSchemeOdeBase(DihuContext context, std::string name) :
  TimeSteppingScheme(context), initialized_(false)
  {
    // get python config
    PythonConfig topLevelSettings = this->context_.getPythonConfig();
    this->specificSettings_ = PythonConfig(topLevelSettings, name);
    
    // initialize output writers
    this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
  }
  
  template<typename FunctionSpaceType, int nComponents>
  Data::TimeStepping<FunctionSpaceType, nComponents> &TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
  data()
  {
    return *data_;
  }
  
  template<typename FunctionSpaceType, int nComponents>
  void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
  setInitialValues()
  {
    // set initial values as given in settings, or set to zero if not given
    std::vector<double> localValues;
    
    bool inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
    if (inputMeshIsGlobal)
    {
      assert(this->data_);
      assert(this->data_->functionSpace());
      const int nDofsGlobal = this->data_->functionSpace()->nDofsGlobal();
      LOG(DEBUG) << "setInitialValues, nDofsGlobal = " << nDofsGlobal;
      
      this->specificSettings_.getOptionVector("initialValues", nDofsGlobal, localValues);
      
      this->data_->functionSpace()->meshPartition()->extractLocalDofsWithoutGhosts(localValues);
    }
    else 
    {
      const int nDofsLocal = this->data_->functionSpace()->nDofsLocalWithoutGhosts();
      this->specificSettings_.getOptionVector("initialValues", nDofsLocal, localValues);
    }
    VLOG(1) << "set initial values to " << localValues;
    
    // set the first component of the solution variable by the given values
    data_->solution()->setValuesWithoutGhosts(0, localValues);
    
    VLOG(1) << data_->solution();
  }
  /*
   * t *emplate<typename DiscretizableInTimeType>
   * std::shared_ptr<typename TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::Data::FieldVariableType> TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
   * solution()
   * {
   * return data_->solution();
}*/
  
  template<typename FunctionSpaceType, int nComponents>
  typename TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::TransferableSolutionDataType TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
  getSolutionForTransferInOperatorSplitting()
  {
    return data_->getSolutionForTransferInOperatorSplitting();
  }
  
  template<typename FunctionSpaceType, int nComponents>
  void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
  setRankSubset(Partition::RankSubset rankSubset)
  {
    data_->setRankSubset(rankSubset);
    //discretizableInTime_.setRankSubset(rankSubset);
  } 
  
  template<typename FunctionSpaceType, int nComponents>
  void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
  reset()
  {
    TimeSteppingScheme::reset();
    //discretizableInTime_.reset();
    
    initialized_ = false;
  }
  
  template<typename FunctionSpaceType, int nComponents>
  void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
  initialize()
  {
    if (initialized_)
      return;
    
    TimeSteppingScheme::initialize();
    LOG(TRACE) << "TimeSteppingSchemeOdeBase::initialize";
    
    //std::shared_ptr<typename DiscretizableInTimeType::FunctionSpace> functionSpace
    //= discretizableInTime_.functionSpace();
    
    //assert(functionSpace->meshPartition());   // assert that the function space was retrieved correctly
    //data_->setFunctionSpace(functionSpace);
    
    // create the vectors in the data object
    //data_->initialize();
    
    // set initial values from settings
    
    // if it did not initialize it,
    // load initial values from config under the timestepping section
    //this->setInitialValues();
    
    //data_->print();
    
    //initialized_ = true;
  }
  
  template<typename FunctionSpaceType, int nComponents>
  void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
  run()
  {
    // initialize
    this->initialize();
    
    // do simulations
    this->advanceTimeSpan();
  }
  
///  template<typename FunctionSpaceType, int nComponents>
///  bool TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
///  knowsMeshType()
///  {
///    return false;
///  }
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///template<typename DiscretizableInTimeType>
///TimeSteppingSchemeOdeBase<DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>::
///TimeSteppingSchemeOdeBase(DihuContext context, std::string name) :
///  TimeSteppingScheme(context), discretizableInTime_(context[name]), initialized_(false)
///{
  // get python config
  ///PythonConfig topLevelSettings = this->context_.getPythonConfig();
  ///this->specificSettings_ = PythonConfig(topLevelSettings, name);

  // initialize output writers
  //this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  // create dirichlet Boundary conditions object
  ///this->dirichletBoundaryConditions_ = std::make_shared<
    ///SpatialDiscretization::DirichletBoundaryConditions<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>
  ///>();
///}
/*
template<typename DiscretizableInTimeType>
Data::TimeStepping<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()> &TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
data()
{
  return *data_;
}
*/
/*
template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
setInitialValues()
{
  // set initial values as given in settings, or set to zero if not given
  std::vector<double> localValues;
  
  bool inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    assert(this->data_);
    assert(this->data_->functionSpace());
    const int nDofsGlobal = this->data_->functionSpace()->nDofsGlobal();
    LOG(DEBUG) << "setInitialValues, nDofsGlobal = " << nDofsGlobal;

    this->specificSettings_.getOptionVector("initialValues", nDofsGlobal, localValues);

    this->data_->functionSpace()->meshPartition()->extractLocalDofsWithoutGhosts(localValues);
  }
  else 
  {
    const int nDofsLocal = this->data_->functionSpace()->nDofsLocalWithoutGhosts();
    this->specificSettings_.getOptionVector("initialValues", nDofsLocal, localValues);
  }
  VLOG(1) << "set initial values to " << localValues;

  // set the first component of the solution variable by the given values
  data_->solution()->setValuesWithoutGhosts(0, localValues);

  VLOG(1) << data_->solution();
}
*/
/*
template<typename DiscretizableInTimeType>
std::shared_ptr<typename TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::Data::FieldVariableType> TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
solution()
{
  return data_->solution();
}*/
/*
template<typename DiscretizableInTimeType>
typename TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::TransferableSolutionDataType TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
getSolutionForTransferInOperatorSplitting()
{
  return data_->getSolutionForTransferInOperatorSplitting();
}
*/
template<typename DiscretizableInTimeType>
DiscretizableInTimeType &TimeSteppingSchemeOde<DiscretizableInTimeType>::
discretizableInTime()
{
  return this->discretizableInTime_;
}
/*
template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
setRankSubset(Partition::RankSubset rankSubset)
{
  data_->setRankSubset(rankSubset);
  discretizableInTime_.setRankSubset(rankSubset);
} 
*/

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOde<DiscretizableInTimeType>::
reset()
{
  TimeSteppingScheme::reset();
  discretizableInTime_.reset();
  
  initialized_ = false;
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOde<DiscretizableInTimeType>::
initialize()
{
  if (initialized_)
    return;
 
  TimeSteppingScheme::initialize();
  LOG(TRACE) << "TimeSteppingSchemeOdeBase::initialize";

  // disable boundary condition handling in finite element method, because Dirichlet BC have to be handled in the system matrix here
  discretizableInTime_.setBoundaryConditionHandlingEnabled(false);

  // initialize underlying DiscretizableInTime object, also with time step width
  discretizableInTime_.initialize();
  discretizableInTime_.initializeForImplicitTimeStepping();   // this performs extra initialization for implicit timestepping methods, i.e. it sets the inverse lumped mass matrix

  std::shared_ptr<typename DiscretizableInTimeType::FunctionSpace> functionSpace
    = discretizableInTime_.functionSpace();

  assert(functionSpace->meshPartition());   // assert that the function space was retrieved correctly
  this->data_->setFunctionSpace(functionSpace);
  
  // set component names in data
  std::vector<std::string> componentNames;
  discretizableInTime_.getComponentNames(componentNames);
  this->data_->setComponentNames(componentNames);
  
  // create the vectors in the data object
  this->data_->initialize();

  // parse boundary conditions, needs functionSpace set
  // initialize dirichlet boundary conditions object which parses dirichlet boundary condition dofs and values from config
  this->dirichletBoundaryConditions_->initialize(this->specificSettings_, this->data_->functionSpace());

  // set initial values from settings

  // load initial values as specified in config under the "CellML" section
  if (!discretizableInTime_.setInitialValues(this->data_->solution()))
  {
    LOG(DEBUG) << "initial values were not set by DiscretizableInTime, set now";

    // if it did not initialize it,
    // load initial values from config under the timestepping section
    this->setInitialValues();
  }
  else
  {
    LOG(DEBUG) << "initial values were set by DiscretizableInTime";
  }
  
  this->data_->print();
  
  initialized_ = true;
}

template<int nStates, typename FunctionSpaceType>
void TimeSteppingSchemeOde<CellmlAdapter<nStates, FunctionSpaceType>>::
initialize()
{
  TimeSteppingSchemeOdeBase<FunctionSpaceType,nStates>::initialize();
  double prefactor = this->discretizableInTime_.prefactor();
  int outputComponentNo = this->discretizableInTime_.outputStateIndex();

  LOG(DEBUG) << "set CellML prefactor=" << prefactor << ", outputComponentNo=" << outputComponentNo;

  this->data_->setPrefactor(prefactor);
  this->data_->setOutputComponentNo(outputComponentNo);
}
/*
template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
run()
{
  // initialize
  this->initialize();

  // do simulations
  this->advanceTimeSpan();
}
*/
template<typename DiscretizableInTimeType>
bool TimeSteppingSchemeOde<DiscretizableInTimeType>::
knowsMeshType()
{
  return this->discretizableInTime_.knowsMeshType();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////


} // namespace
