#include "time_stepping_scheme/time_stepping_scheme_ode.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::TimeSteppingSchemeOdeBaseDiscretizable(DihuContext context, std::string name) :
  TimeSteppingSchemeOdeBase<typename DiscretizableInTimeType::FunctionSpace,DiscretizableInTimeType::nComponents()>::
  TimeSteppingSchemeOdeBase(context, name), discretizableInTime_(this->context_), initialized_(false)
{
  //initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  //create dirichlet Boundary conditions object
  this->dirichletBoundaryConditions_ = std::make_shared<
    SpatialDiscretization::DirichletBoundaryConditions<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>
  >(this->context_);
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::
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
  this->data_->solution()->setValuesWithoutGhosts(0, localValues);

  VLOG(1) << this->data_->solution();
}

template<typename DiscretizableInTimeType>
DiscretizableInTimeType &TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::
discretizableInTime()
{
  return this->discretizableInTime_;
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::
setRankSubset(Partition::RankSubset rankSubset)
{
  TimeSteppingSchemeOdeBase<typename DiscretizableInTimeType::Functionspace,DiscretizableInTimeType::nComponents()>::
  setRankSubset(rankSubset);
  discretizableInTime_.setRankSubset(rankSubset);
} 

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::
reset()
{
  TimeSteppingSchemeOdeBase<typename DiscretizableInTimeType::FunctionSpace,DiscretizableInTimeType::nComponents()>::
  reset();
  discretizableInTime_.reset();
  this->data_.reset();
  
  initialized_ = false;
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::
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

  // retrieve the function space from the discretizable in time object, this is used for the data object
  std::shared_ptr<typename DiscretizableInTimeType::FunctionSpace> functionSpace
    = discretizableInTime_.functionSpace();

  assert(functionSpace->meshPartition());   // assert that the function space was retrieved correctly
  assert(this->data_);
  this->data_->setFunctionSpace(functionSpace);
  
  // set component names in data
  std::vector<std::string> componentNames;
  discretizableInTime_.getComponentNames(componentNames);
  this->data_->setComponentNames(componentNames);
  
  // create the vectors in the data object
  this->data_->initialize();

  // parse boundary conditions, needs functionSpace set
  // initialize dirichlet boundary conditions object which parses dirichlet boundary condition dofs and values from config
  this->dirichletBoundaryConditions_->initialize(this->specificSettings_, this->data_->functionSpace(), "dirichletBoundaryConditions");
  //TODO: add Neumann BC

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
  VLOG(1) << "initial solution vector: " << *this->data_->solution();
  
  //output initial values
  this->outputWriterManager_.writeOutput(*this->data_, 0, 0);
  
  this->data_->print();
  
  initialized_ = true;
}

template<typename DiscretizableInTimeType>
std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<typename DiscretizableInTimeType::FunctionSpace,DiscretizableInTimeType::nComponents()>>
TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::
dirichletBoundaryConditions()
{
  return dirichletBoundaryConditions_;
}

template<int nStates, typename FunctionSpaceType>
void TimeSteppingSchemeOde<CellmlAdapter<nStates, FunctionSpaceType>>::
initialize()
{
  TimeSteppingSchemeOdeBaseDiscretizable<CellmlAdapter<nStates, FunctionSpaceType>>::initialize();
  double prefactor = this->discretizableInTime_.prefactor();
  int outputComponentNo = this->discretizableInTime_.outputStateIndex();

  LOG(DEBUG) << "set CellML prefactor=" << prefactor << ", outputComponentNo=" << outputComponentNo;

  this->data_->setPrefactor(prefactor);
  this->data_->setOutputComponentNo(outputComponentNo);
}

template<typename DiscretizableInTimeType>
bool TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::
knowsMeshType()
{
  return this->discretizableInTime_.knowsMeshType();
}

} // namespace
