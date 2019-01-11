#include "time_stepping_scheme/time_stepping_scheme_ode.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
TimeSteppingSchemeOdeBase(DihuContext context, std::string name) :
  TimeSteppingScheme(context), discretizableInTime_(context[name]), initialized_(false)
{
  // get python config
  PythonConfig topLevelSettings = this->context_.getPythonConfig();
  this->specificSettings_ = PythonConfig(topLevelSettings, name);

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  // create dirichlet Boundary conditions object
  this->dirichletBoundaryConditions_ = std::make_shared<
    SpatialDiscretization::DirichletBoundaryConditions<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>
  >();
}

template<typename DiscretizableInTimeType>
Data::TimeStepping<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()> &TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
data()
{
  return *data_;
}

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

template<typename DiscretizableInTimeType>
typename TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::TransferableSolutionDataType TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
getSolutionForTransfer()
{
  return data_->getSolutionForTransfer();
}

template<typename DiscretizableInTimeType>
DiscretizableInTimeType &TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
discretizableInTime()
{
  return this->discretizableInTime_;
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
setRankSubset(Partition::RankSubset rankSubset)
{
  data_->setRankSubset(rankSubset);
  discretizableInTime_.setRankSubset(rankSubset);
} 
 
template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
reset()
{
  TimeSteppingScheme::reset();
  discretizableInTime_.reset();
  
  initialized_ = false;
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
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
  data_->setFunctionSpace(functionSpace);
  
  // set component names in data
  std::vector<std::string> componentNames;
  discretizableInTime_.getComponentNames(componentNames);
  data_->setComponentNames(componentNames);
  
  // create the vectors in the data object
  data_->initialize();

  // parse boundary conditions, needs functionSpace set
  // initialize dirichlet boundary conditions object which parses dirichlet boundary condition dofs and values from config
  this->dirichletBoundaryConditions_->initialize(this->specificSettings_, this->data_->functionSpace());

  // set initial values from settings

  // load initial values as specified in config under the "CellML" section
  if (!discretizableInTime_.setInitialValues(data_->solution()))
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
  
  data_->print();
  
  initialized_ = true;
}

template<int nStates, typename FunctionSpaceType>
void TimeSteppingSchemeOde<CellmlAdapter<nStates, FunctionSpaceType>>::
initialize()
{
  TimeSteppingSchemeOdeBase<CellmlAdapter<nStates, FunctionSpaceType>>::initialize();
  double prefactor = this->discretizableInTime_.prefactor();
  int outputComponentNo = this->discretizableInTime_.outputStateIndex();

  LOG(DEBUG) << "set CellML prefactor=" << prefactor << ", outputComponentNo=" << outputComponentNo;

  this->data_->setPrefactor(prefactor);
  this->data_->setOutputComponentNo(outputComponentNo);
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
run()
{
  // initialize
  this->initialize();

  // do simulations
  this->advanceTimeSpan();
}

template<typename DiscretizableInTimeType>
bool TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::
knowsMeshType()
{
  return this->discretizableInTime_.knowsMeshType();
}


} // namespace
