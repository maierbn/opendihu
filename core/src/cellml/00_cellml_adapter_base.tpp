#include "cellml/00_cellml_adapter_base.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "data_management/output_connector_data.h"
#include "mesh/mesh_manager/mesh_manager.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
std::array<double,nStates_> CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::statesInitialValues_;

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
bool CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::statesInitialValuesinitialized_ = false;

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
CellmlAdapterBase(DihuContext context) :
  context_(context), specificSettings_(PythonConfig(context_.getPythonConfig(), "CellML")), data_(context_)
{
  outputWriterManager_.initialize(this->context_, specificSettings_);
  LOG(TRACE) << "CellmlAdapterBase constructor";
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
CellmlAdapterBase(DihuContext context, bool noNewOutputWriter) :
  context_(context), specificSettings_(PythonConfig(context_.getPythonConfig(), "CellML"))
{
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
~CellmlAdapterBase()
{
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
constexpr int CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
nComponents()
{
  return nStates_;
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
constexpr int CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
nStates()
{
  return nStates_;
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
setSolutionVariable(std::shared_ptr<FieldVariableStates> states)
{
  this->data_.setStatesVariable(states);
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
setOutputConnectorData(std::shared_ptr<::Data::OutputConnectorData<FunctionSpaceType,nStates_>> outputConnectorDataTimeStepping)
{
  // add all intermediate values for transfer (option "intermediatesForTransfer"), which are stored in this->data_.getOutputConnectorData()
  // at the end of outputConnectorDataTimeStepping

  // loop over intermediates that should be transferred
  for (typename std::vector<::Data::ComponentOfFieldVariable<FunctionSpaceType,nIntermediates_>>::iterator iter
    = this->data_.getOutputConnectorData()->variable2.begin(); iter != this->data_.getOutputConnectorData()->variable2.end(); iter++)
  {
    int componentNo = iter->componentNo;
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nIntermediates_>> values = iter->values;

    values->setRepresentationGlobal();

    // The intermediate field variables have 'nIntermediates_' components, but the field variables in the outputConnectorDataTimeStepping object
    // have only 1 component. Therefore, we create new field variables with 1 components each that reuse the Petsc Vec's of the intermediate field variables.

    // get the parameters to create the new field variable
    std::string name = values->componentName(componentNo);
    const std::vector<std::string> componentNames{values->componentName(componentNo)};
    const bool reuseData = true;

    // create the new field variable with only the one component
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> newFieldVariable
      = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,1>>(*values, name, componentNames, reuseData);

    LOG(DEBUG) << "CellmlAdapterBase::setOutputConnectorData add FieldVariable " << newFieldVariable << " for intermediate " << componentNo << "," << name;

    // add this component to outputConnector of data time stepping
    outputConnectorDataTimeStepping->addFieldVariable2(newFieldVariable);
  }
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
initialize()
{
  LOG(TRACE) << "CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::initialize";

  if (VLOG_IS_ON(1))
  {
    LOG(DEBUG) << "CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::initialize querying meshManager for mesh";
    LOG(DEBUG) << "specificSettings_: ";
    PythonUtility::printDict(specificSettings_.pyObject());
  }
  
  // add this solver to the solvers diagram, which is a SVG file that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("CellmlAdapter");

  // create a mesh if there is not yet one assigned, function space FunctionSpace::Generic
  if (!functionSpace_)
  {
    functionSpace_ = context_.meshManager()->functionSpace<FunctionSpaceType>(specificSettings_);  // create initialized mesh
  }
  LOG(DEBUG) << "Cellml mesh has " << functionSpace_->nNodesLocalWithoutGhosts() << " local nodes";

  //store number of instances
  nInstances_ = functionSpace_->nNodesLocalWithoutGhosts();


  // initialize source code generator
  std::string modelFilename = this->specificSettings_.getOptionString("modelFilename", "");
  if (this->specificSettings_.hasKey("sourceFilename"))
  {
    LOG(WARNING) << this->specificSettings_ << " Option \"sourceFilename\" has been renamed to \"modelFilename\".";
  }

  // add explicitely defined parameters that replace intermediates and constants
  std::vector<int> parametersUsedAsIntermediate;  ///< explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array
  std::vector<int> parametersUsedAsConstant;  ///< explicitely defined parameters that will be copied to constants, this vector contains the indices of the constants

  this->specificSettings_.getOptionVector("parametersUsedAsIntermediate", parametersUsedAsIntermediate);
  this->specificSettings_.getOptionVector("parametersUsedAsConstant", parametersUsedAsConstant);

  // initialize parameters
  std::vector<double> parametersInitialValues;
  if (this->specificSettings_.hasKey("parametersInitialValues"))
  {
    LOG(DEBUG) << "load parametersInitialValues from config";
    specificSettings_.getOptionVector("parametersInitialValues", parametersInitialValues);
  }
  else
  {
    LOG(DEBUG) << "Config does not contain key \"parametersInitialValues\"";
  }

  cellmlSourceCodeGenerator_.initialize(modelFilename, nInstances_, nStates_, nIntermediates_,
                                        parametersUsedAsIntermediate, parametersUsedAsConstant, parametersInitialValues);

  // initialize data, i.e. states and intermediates field variables
  data_.setFunctionSpace(functionSpace_);
  data_.setIntermediateNames(cellmlSourceCodeGenerator_.intermediateNames());
  data_.initialize();

  initializeStatesToEquilibrium_ = this->specificSettings_.getOptionBool("initializeStatesToEquilibrium", false);
  if (initializeStatesToEquilibrium_)
  {
    initializeStatesToEquilibriumTimestepWidth_ = this->specificSettings_.getOptionDouble("initializeStatesToEquilibriumTimestepWidth", 1e-4);
  }
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
template<typename FunctionSpaceType2>
bool CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nStates_>> initialValues)
{
  LOG(TRACE) << "CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::setInitialValues";

  if (!statesInitialValuesinitialized_)
  {

    // initialize states
    if (this->specificSettings_.hasKey("statesInitialValues"))
    {
      LOG(DEBUG) << "set initial values from config";

      // statesInitialValues gives the initial state values for one instance of the problem. it is used for all instances.
      std::array<double,nStates_> statesInitialValuesFromConfig = this->specificSettings_.template getOptionArray<double,nStates_>("statesInitialValues", 0);

      // store initial values to statesInitialValues_
      std::copy(statesInitialValuesFromConfig.begin(), statesInitialValuesFromConfig.end(), statesInitialValues_.begin());
    }
    else
    {
      LOG(DEBUG) << "set initial values from source file";

      // parsing the source file was already done
      // get initial values from source code generator
      std::vector<double> statesInitialValuesGenerator = cellmlSourceCodeGenerator_.statesInitialValues();
      assert(statesInitialValuesGenerator.size() == nStates_);

      std::copy(statesInitialValuesGenerator.begin(), statesInitialValuesGenerator.end(), statesInitialValues_.begin());
    }

    if (initializeStatesToEquilibrium_)
    {
      this->initializeToEquilibriumValues(statesInitialValues_);
    }
    statesInitialValuesinitialized_ = true;
  }

  // Here we have the initial values for the states in the statesInitialValues vector, only for one instance.
  VLOG(1) << "statesInitialValues: " << statesInitialValues_;
  const std::vector<std::array<double,nStates_>> statesAllInstances(nInstances_, statesInitialValues_);

  VLOG(1) << "statesAllInstances: " << statesAllInstances << ", nInstances: " << nInstances_ << ", nStates_ per instances: " << statesInitialValues_.size();
  initialValues->setValuesWithoutGhosts(statesAllInstances);

  VLOG(1) << "initialValues: " << *initialValues;
  return true;
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
functionSpace()
{
  return functionSpace_;
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
getNumbers(int& nInstances, int& nIntermediates, int& nParameters)
{
  nInstances = nInstances_;
  nIntermediates = nIntermediates_;
  nParameters = cellmlSourceCodeGenerator_.nParameters();
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
getStatesIntermediatesForTransfer(std::vector<int> &statesForTransfer, std::vector<int> &intermediatesForTransfer)
{
  data_.getStatesIntermediatesForTransfer(statesForTransfer, intermediatesForTransfer);
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
getStateNames(std::vector<std::string> &stateNames)
{
  stateNames = this->cellmlSourceCodeGenerator_.stateNames();
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
typename CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::Data &CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
data()
{
  return this->data_;
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
CellmlSourceCodeGenerator &CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
cellmlSourceCodeGenerator()
{
  return this->cellmlSourceCodeGenerator_;
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
PythonConfig CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
specificSettings()
{
  return this->specificSettings_;
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
constexpr int  CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
nIntermediates() const
{
  return nIntermediates_;
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
std::shared_ptr<typename CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::OutputConnectorDataType>
CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
getOutputConnectorData()
{
  return this->data_.getOutputConnectorData();
}
