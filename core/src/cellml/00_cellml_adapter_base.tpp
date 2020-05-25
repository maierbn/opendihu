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
  context_(context), specificSettings_(PythonConfig(context_.getPythonConfig(), "CellML")),
  data_(context_), cellmlSourceCodeGenerator_()
{
  outputWriterManager_.initialize(this->context_, specificSettings_);
  LOG(TRACE) << "CellmlAdapterBase constructor";
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
CellmlAdapterBase(DihuContext context, bool initializeOutputWriter) :
  context_(context), specificSettings_(PythonConfig(context_.getPythonConfig(), "CellML")),
  data_(context_), cellmlSourceCodeGenerator_()
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
  // add all state and intermediate values for transfer (option "intermediatesForTransfer"), which are stored in this->data_.getOutputConnectorData()
  // at the end of outputConnectorDataTimeStepping

  // loop over states that should be transferred
  for (typename std::vector<::Data::ComponentOfFieldVariable<FunctionSpaceType,nStates_>>::iterator iter
    = this->data_.getOutputConnectorData()->variable1.begin(); iter != this->data_.getOutputConnectorData()->variable1.end(); iter++)
  {
    // skip the first "states" entry of statesToTransfer, because this is the solution variable of the timestepping scheme and therefore
    // the timestepping scheme has already added it to the outputConnectorDataTimeStepping object
    if (iter == this->data_.getOutputConnectorData()->variable1.begin())
      continue;

    int componentNo = iter->componentNo;
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStates_>> values = iter->values;

    values->setRepresentationGlobal();

    // The state field variables have 'nStates_' components and can be reused.
    std::string name = values->componentName(componentNo);
    LOG(DEBUG) << "CellmlAdapterBase::setOutputConnectorData add FieldVariable " << *values << " for state " << componentNo << "," << name;

    // add this component to outputConnector of data time stepping
    outputConnectorDataTimeStepping->addFieldVariable(values, componentNo);
  }

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

  // initialize parameters
  std::vector<double> parametersInitialValues;
  specificSettings_.getOptionVector("parametersInitialValues", parametersInitialValues);

  // parse mappings, which initializes the following variables:
  // parametersUsedAsIntermediate, parametersUsedAsConstant,

  // add explicitely defined parameters that replace intermediates and constants
  std::vector<int> parametersUsedAsIntermediate;  //< explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array
  std::vector<int> parametersUsedAsConstant;      //< explicitely defined parameters that will be copied to constants, this vector contains the indices of the constants

  std::vector<int> &statesForTransfer = data_.statesForTransfer();
  std::vector<int> &intermediatesForTransfer = data_.intermediatesForTransfer();
  std::vector<int> &parametersForTransfer = data_.parametersForTransfer();

  // parse the source code and initialize the names of states, intermediates and constants, which are needed for initializeMappings
  cellmlSourceCodeGenerator_.initializeNames(modelFilename, nInstances_, nStates_, nIntermediates_);

  // initialize all information from python settings key "mappings", this sets parametersUsedAsIntermediates/States and outputIntermediate/StatesIndex
  initializeMappings(parametersUsedAsIntermediate, parametersUsedAsConstant,
                     statesForTransfer, intermediatesForTransfer, parametersForTransfer);

  // initialize data, i.e. states and intermediates field variables
  data_.setFunctionSpace(functionSpace_);
  data_.setIntermediateNames(cellmlSourceCodeGenerator_.intermediateNames());
  data_.initialize();

  // get the data_.parameters() raw pointer
  data_.prepareParameterValues();

  // parse the source code completely and store source code, needs data initialized in order to store initial parameter values
  cellmlSourceCodeGenerator_.initializeSourceCode(parametersUsedAsIntermediate, parametersUsedAsConstant, parametersInitialValues, nIntermediates_, data_.parameterValues());

  // restore the raw pointer of data_.parameters()
  data_.restoreParameterValues();

  initializeStatesToEquilibrium_ = this->specificSettings_.getOptionBool("initializeStatesToEquilibrium", false);
  if (initializeStatesToEquilibrium_)
  {
    initializeStatesToEquilibriumTimestepWidth_ = this->specificSettings_.getOptionDouble("initializeStatesToEquilibriumTimestepWidth", 1e-4);
  }
}


template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
initializeMappings(std::vector<int> &parametersUsedAsIntermediate, std::vector<int> &parametersUsedAsConstant,
                   std::vector<int> &statesForTransfer, std::vector<int> &intermediatesForTransfer, std::vector<int> &parametersForTransfer)
{
  if (this->specificSettings_.hasKey("mappings"))
  {
    // example for mappings:
    /*
    "mappings": {
      ("parameter",0): ("state",0),
      ("parameter",1): ("intermediate",0),
      ("parameter",2): ("constant",0),
      ("outputConnectorSlot",0): ("parameter",0),
    }
    */
    const std::vector<std::string> &intermediateNames = cellmlSourceCodeGenerator_.intermediateNames();
    const std::vector<std::string> &stateNames = cellmlSourceCodeGenerator_.stateNames();
    const std::vector<std::string> &constantNames = cellmlSourceCodeGenerator_.constantNames();

    using EntryTypePython = std::pair<std::string,PyObject *>;

    std::pair<EntryTypePython,EntryTypePython> item = this->specificSettings_.template getOptionDictBegin<EntryTypePython,EntryTypePython>("mappings");

    for (;!this->specificSettings_.getOptionDictEnd("mappings");
         this->specificSettings_.template getOptionDictNext<EntryTypePython,EntryTypePython>("mappings",item))
    {
      EntryTypePython entries[2] = {item.first, item.second};
      std::vector<std::pair<std::string,int>> entriesWithNos;

      // parse the strings in the entries and convert them to stateNo, intermediateNo or constantNo
      for (EntryTypePython entry : entries)
      {
        if (entry.first == "parameter" || entry.first == "outputConnectorSlot")
        {
          int parameterNo = PythonUtility::convertFromPython<int>::get(entry.second);
          entriesWithNos.push_back(std::make_pair(entry.first, parameterNo));
        }
        else if (entry.first == "state")
        {
          int stateNo = 0;
          std::string stateName = PythonUtility::convertFromPython<std::string>::get(entry.second);

          // if parsed string is a number, this is directly the state no
          if (!stateName.empty() && stateName.find_first_not_of("0123456789") == std::string::npos)
          {
            stateNo = atoi(stateName.c_str());
          }
          else
          {
            // find if the string is a state name given in stateNames
            const std::vector<std::string>::const_iterator pos = std::find(stateNames.begin(), stateNames.end(), stateName);
            if (pos != stateNames.end())
            {
              stateNo = std::distance(stateNames.begin(),pos);
            }
            else
            {
              LOG(ERROR) << "In " << this->specificSettings_ << "[\"mappings\"], state \"" << stateName << "\" is not valid. "
                << "Valid state names are: " << stateNames << ".\n"
                << "You can also specify the state no. instead of the name.";
            }
          }

          entriesWithNos.push_back(std::make_pair(entry.first, stateNo));
        }
        else if (entry.first == "intermediate")
        {
          int intermediateNo = 0;
          std::string intermediateName = PythonUtility::convertFromPython<std::string>::get(entry.second);

          // if parsed string is a number, this is directly the intermediate no
          if (!intermediateName.empty() && intermediateName.find_first_not_of("0123456789") == std::string::npos)
          {
            intermediateNo = atoi(intermediateName.c_str());
          }
          else
          {
            // find if the string is a intermediate name given in intermediateNames
            std::vector<std::string>::const_iterator pos = std::find(intermediateNames.begin(), intermediateNames.end(), intermediateName);
            if (pos != intermediateNames.end())
            {
              intermediateNo = std::distance(intermediateNames.begin(),pos);
            }
            else
            {
              LOG(ERROR) << "In " << this->specificSettings_ << "[\"mappings\"], intermediate \"" << intermediateName << "\" is not valid. "
                << "Valid intermediate names are: " << intermediateNames << ".\n"
                << "You can also specify the intermediate no. instead of the name.";
            }
          }

          entriesWithNos.push_back(std::make_pair(entry.first, intermediateNo));
        }
        else if (entry.first == "constant")
        {
          int constantNo = 0;
          std::string constantName = PythonUtility::convertFromPython<std::string>::get(entry.second);

          // if parsed string is a number, this is directly the constant no
          if (!constantName.empty() && constantName.find_first_not_of("0123456789") == std::string::npos)
          {
            constantNo = atoi(constantName.c_str());
          }
          else
          {
            // find if the string is a constant name given in constantNames
            std::vector<std::string>::const_iterator pos = std::find(constantNames.begin(), constantNames.end(), constantName);
            if (pos != constantNames.end())
            {
              constantNo = std::distance(constantNames.begin(),pos);
            }
            else
            {
              LOG(ERROR) << "In " << this->specificSettings_ << "[\"mappings\"], constant \"" << constantName << "\" is not valid. "
                << "Valid constant names are: " << constantNames << ".\n"
                << "You can also specify the constant no. instead of the name.";
            }
          }

          entriesWithNos.push_back(std::make_pair(entry.first, constantNo));
        }
        else
        {
          LOG(FATAL) << "In " << this->specificSettings_ << "[\"mappings\"], invalid field \"" << entry.first << "\", allowed values are: "
            << "\"parameter\", \"outputConnectorSlot\", \"state\", \"intermediate\", \"constant\".";
        }
      }

      LOG(DEBUG) << "parsed next item, entriesWithNos: " << entriesWithNos;

      // make sure that the first part is "paramater" or "outputConnectorSlot"
      if (entriesWithNos[1].first == "parameter" || entriesWithNos[1].first == "outputConnectorSlot")
      {
        std::swap(entriesWithNos[0], entriesWithNos[1]);
      }

      // handle parameter
      if (entriesWithNos[0].first == "parameter")
      {
        int parameterNo = entriesWithNos[0].second;    // the no. of the parameter
        int fieldNo = entriesWithNos[1].second;        // the no. of the state, intermediate or constant

        LOG(DEBUG) << "parameter " << parameterNo << " to \"" << entriesWithNos[1].first << "\", fieldNo " << fieldNo;

        if (entriesWithNos[1].first == "intermediate")
        {
          // add no of the constant to parametersUsedAsIntermediate
          if (parametersUsedAsIntermediate.size() <= parameterNo)
          {
            parametersUsedAsIntermediate.resize(parameterNo+1, -1);
          }
          parametersUsedAsIntermediate[parameterNo] = fieldNo;
        }
        else if (entriesWithNos[1].first == "constant")
        {
          // add no of the constant to parametersUsedAsConstant
          if (parametersUsedAsConstant.size() <= parameterNo)
          {
            parametersUsedAsConstant.resize(parameterNo+1, -1);
          }
          parametersUsedAsConstant[parameterNo] = fieldNo;
        }
        else
        {
          LOG(ERROR) << "In " << this->specificSettings_ << "[\"mappings\"], you can map "
            << "paramaters only to intermediates and constants, not states.";
        }
      }
      else
      {
        // output connector slots

        int slotNo = entriesWithNos[0].second;    // the no. of the output connector slot
        int fieldNo = entriesWithNos[1].second;        // the no. of the state, intermediate or constant

        LOG(DEBUG) << "output connector slot " << slotNo << " to \"" << entriesWithNos[1].first << "\", fieldNo " << fieldNo;

        if (entriesWithNos[1].first == "state")
        {
          // add no of the constant to statesForTransfer
          if (statesForTransfer.size() <= slotNo)
          {
            statesForTransfer.resize(slotNo+1, -1);
          }
          statesForTransfer[slotNo] = fieldNo;
        }
        else if (entriesWithNos[1].first == "intermediate")
        {
          // add no of the constant to intermediatesForTransfer
          if (intermediatesForTransfer.size() <= slotNo)
          {
            intermediatesForTransfer.resize(slotNo+1, -1);
          }
          intermediatesForTransfer[slotNo] = fieldNo;
        }
        else if (entriesWithNos[1].first == "parameter")
        {
          // add no of the constant to parametersForTransfer
          if (parametersForTransfer.size() <= slotNo)
          {
            parametersForTransfer.resize(slotNo+1, -1);
          }
          parametersForTransfer[slotNo] = fieldNo;
        }
        else
        {
          LOG(ERROR) << "In " << this->specificSettings_ << "[\"mappings\"], you can map "
            << "outputConnectorSlots only to states, intermediates and parameters, not constants.";
        }
      }
    }

    // reorder slot nos such that states go first, then intermediates, then parameters
    std::vector<int> temporaryVector;
    std::copy_if(statesForTransfer.begin(), statesForTransfer.end(),
                 std::back_insert_iterator<std::vector<int>>(temporaryVector), [](int a){return a != -1;});
    statesForTransfer = temporaryVector;

    temporaryVector.clear();
    std::copy_if(intermediatesForTransfer.begin(), intermediatesForTransfer.end(),
                 std::back_insert_iterator<std::vector<int>>(temporaryVector), [](int a){return a != -1;});
    intermediatesForTransfer = temporaryVector;

    temporaryVector.clear();
    std::copy_if(parametersForTransfer.begin(), parametersForTransfer.end(),
                 std::back_insert_iterator<std::vector<int>>(temporaryVector), [](int a){return a != -1;});
    parametersForTransfer = temporaryVector;

    // create a string of parameter and contstant mappings
    std::stringstream s;
    for (int i = 0; i < parametersUsedAsIntermediate.size(); i++)
    {
      int intermediateNo = parametersUsedAsIntermediate[i];
      if (intermediateNo != -1 && intermediateNo < nIntermediates_)
      {
        s << "  parameter " << i << " maps to intermediate " << intermediateNo << " (\"" << intermediateNames[intermediateNo] << "\")\n";
      }
    }
    for (int i = 0; i < parametersUsedAsConstant.size(); i++)
    {
      int constantNo = parametersUsedAsConstant[i];
      if (constantNo != -1 && constantNo < constantNames.size())
      {
        s << "  parameter " << i << " maps to constant " << constantNo << " (\"" << constantNames[constantNo] << "\")\n";
      }
    }

    // check that parameters are all mapped
    for (int i = 0; i < parametersUsedAsIntermediate.size(); i++)
    {
      int intermediateNo = parametersUsedAsIntermediate[i];
      if (intermediateNo == -1 && intermediateNo < nIntermediates_)
      {
        LOG(ERROR) << "Parameter no. " << i << " is not mapped to any intermediate. "
          << "(Note that parameters have to be ordered in a way that the first are mapped to intermediates and the last to constants.)\n"
          << " Parsed mapping:\n" << s.str();
      }
    }

    for (int i = 0; i < parametersUsedAsConstant.size(); i++)
    {
      int constantNo = parametersUsedAsConstant[i];
      if (constantNo == -1 && constantNo < constantNames.size())
      {
        LOG(ERROR) << "Parameter no. " << parametersUsedAsIntermediate.size() + i << " is not mapped to any constant. "
          << "(Note that parameters have to be ordered in a way that the first are mapped to intermediates and the last to constants.)\n"
          << " Parsed mapping:\n" << s.str();
      }
    }
  }
  else
  {
    this->specificSettings_.getOptionVector("parametersUsedAsIntermediate", parametersUsedAsIntermediate);
    this->specificSettings_.getOptionVector("parametersUsedAsConstant", parametersUsedAsConstant);

    this->specificSettings_.template getOptionVector<int>("statesForTransfer", statesForTransfer);
    this->specificSettings_.template getOptionVector<int>("intermediatesForTransfer", intermediatesForTransfer);
    this->specificSettings_.template getOptionVector<int>("parametersForTransfer", parametersForTransfer);
  }

  LOG(DEBUG) << "parametersUsedAsIntermediate: " << parametersUsedAsIntermediate;
  LOG(DEBUG) << "parametersUsedAsConstant:     " << parametersUsedAsConstant;
  LOG(DEBUG) << "statesForTransfer:            " << statesForTransfer;
  LOG(DEBUG) << "intermediatesForTransfer:     " << intermediatesForTransfer;
  LOG(DEBUG) << "parametersForTransfer:        " << parametersForTransfer;

  if (this->specificSettings_.hasKey("outputIntermediateIndex"))
  {
    LOG(WARNING) << specificSettings_ << "[\"outputIntermediateIndex\"] is no longer a valid option, use \"intermediatesForTransfer\" instead!";
  }

  if (this->specificSettings_.hasKey("outputStateIndex"))
  {
    LOG(WARNING) << specificSettings_ << "[\"outputStateIndex\"] is no longer a valid option, use \"statesForTransfer\" instead!";
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
    if (this->specificSettings_.hasKey("statesInitialValues") && !this->specificSettings_.isEmpty("statesInitialValues"))
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
std::vector<int> &CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
statesForTransfer()
{
  return data_.statesForTransfer();
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
std::vector<int> &CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::
intermediatesForTransfer()
{
  return data_.intermediatesForTransfer();
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
