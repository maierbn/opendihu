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

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
std::array<double,nStates_> CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::statesInitialValues_;

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
bool CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::statesInitialValuesinitialized_ = false;

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
CellmlAdapterBase(DihuContext context) :
  context_(context), specificSettings_(PythonConfig(context_.getPythonConfig(), "CellML")),
  data_(context_), cellmlSourceCodeGenerator_()
{
  outputWriterManager_.initialize(this->context_, specificSettings_);
  LOG(TRACE) << "CellmlAdapterBase constructor";
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
CellmlAdapterBase(DihuContext context, bool initializeOutputWriter) :
  context_(context), specificSettings_(PythonConfig(context_.getPythonConfig(), "CellML")),
  data_(context_), cellmlSourceCodeGenerator_()
{
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
~CellmlAdapterBase()
{
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
constexpr int CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
nComponents()
{
  return nStates_;
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
constexpr int CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
nStates()
{
  return nStates_;
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
setSolutionVariable(std::shared_ptr<FieldVariableStates> states)
{
  this->data_.setStatesVariable(states);
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
setOutputConnectorData(std::shared_ptr<::Data::OutputConnectorData<FunctionSpaceType,nStates_>> outputConnectorDataTimeStepping)
{
  // add all state and algebraic values for transfer (option "algebraicsForTransfer"), which are stored in this->data_.getOutputConnectorData()
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
    LOG(DEBUG) << "CellmlAdapterBase::setOutputConnectorData add FieldVariable " << *values << " (" << values->name() << ") for state " << componentNo << "," << name;

    // add this component to outputConnector of data time stepping
    outputConnectorDataTimeStepping->addFieldVariable(values, componentNo);
  }

  // loop over algebraics that should be transferred
  for (typename std::vector<::Data::ComponentOfFieldVariable<FunctionSpaceType,nAlgebraics_>>::iterator iter
    = this->data_.getOutputConnectorData()->variable2.begin(); iter != this->data_.getOutputConnectorData()->variable2.end(); iter++)
  {
    int componentNo = iter->componentNo;
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nAlgebraics_>> values = iter->values;

    values->setRepresentationGlobal();

    // The algebraic field variables have 'nAlgebraics_' components, but the field variables in the outputConnectorDataTimeStepping object
    // have only 1 component. Therefore, we create new field variables with 1 components each that reuse the Petsc Vec's of the algebraic field variables.

    // get the parameters to create the new field variable
    std::string name = values->componentName(componentNo);
    const std::vector<std::string> componentNames{"0"};
    const bool reuseData = true;

    // create the new field variable with only the one component, the component given by componentNo
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> newFieldVariable
      = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,1>>(*values, name, componentNames, reuseData, componentNo);

    LOG(DEBUG) << "CellmlAdapterBase::setOutputConnectorData add FieldVariable2 " << newFieldVariable << " with name " << name << " for component no " << componentNo
      << ", this reuses the data from \"" << values->name() << "\".";

    // add this component to outputConnector of data time stepping
    outputConnectorDataTimeStepping->addFieldVariable2(newFieldVariable);
  }
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
initialize()
{
  LOG(TRACE) << "CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::initialize";

  if (VLOG_IS_ON(1))
  {
    LOG(DEBUG) << "CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::initialize querying meshManager for mesh";
    LOG(DEBUG) << "specificSettings_: ";
    PythonUtility::printDict(specificSettings_.pyObject());
  }
  
  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
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
  // parametersUsedAsAlgebraic, parametersUsedAsConstant,

  // add explicitely defined parameters that replace algebraics and constants
  std::vector<int> parametersUsedAsAlgebraic;  //< explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array
  std::vector<int> parametersUsedAsConstant;      //< explicitely defined parameters that will be copied to constants, this vector contains the indices of the constants

  std::vector<int> &statesForTransfer = data_.statesForTransfer();
  std::vector<int> &algebraicsForTransfer = data_.algebraicsForTransfer();
  std::vector<int> &parametersForTransfer = data_.parametersForTransfer();

  std::vector<std::string> parameterNames;

  // parse the source code and initialize the names of states, algebraics and constants, which are needed for initializeMappings
  cellmlSourceCodeGenerator_.initializeNames(modelFilename, nInstances_, nStates_, nAlgebraics_);

  // initialize all information from python settings key "mappings", this sets parametersUsedAsAlgebraics/States and outputAlgebraic/StatesIndex
  initializeMappings(parametersUsedAsAlgebraic, parametersUsedAsConstant,
                     statesForTransfer, algebraicsForTransfer, parametersForTransfer, parameterNames);

  // initialize data, i.e. states and algebraics field variables
  data_.setFunctionSpace(functionSpace_);
  data_.setAlgebraicAndParameterNames(cellmlSourceCodeGenerator_.algebraicNames(), parameterNames);
  data_.initialize();

  // get the data_.parameters() raw pointer
  data_.prepareParameterValues();

  // parse the source code completely and store source code, needs data initialized in order to store initial parameter values
  cellmlSourceCodeGenerator_.initializeSourceCode(parametersUsedAsAlgebraic, parametersUsedAsConstant, parametersInitialValues, nAlgebraics_, data_.parameterValues());

  // restore the raw pointer of data_.parameters()
  data_.restoreParameterValues();

  initializeStatesToEquilibrium_ = this->specificSettings_.getOptionBool("initializeStatesToEquilibrium", false);
  if (initializeStatesToEquilibrium_)
  {
    initializeStatesToEquilibriumTimestepWidth_ = this->specificSettings_.getOptionDouble("initializeStatesToEquilibriumTimestepWidth", 1e-4);
  }
}


template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
initializeMappings(std::vector<int> &parametersUsedAsAlgebraic, std::vector<int> &parametersUsedAsConstant,
                   std::vector<int> &statesForTransfer, std::vector<int> &algebraicsForTransfer, std::vector<int> &parametersForTransfer, std::vector<std::string> &parameterNames)
{
  if (this->specificSettings_.hasKey("mappings"))
  {
    // example for mappings:
    /*
    "mappings": {
      ("parameter",0): ("state",0),
      ("parameter",1): ("algebraic",0),
      ("parameter",2): ("constant",0),
      ("outputConnectorSlot",0): ("parameter",0),
    }
    */
    const std::vector<std::string> &algebraicNames = cellmlSourceCodeGenerator_.algebraicNames();
    const std::vector<std::string> &stateNames = cellmlSourceCodeGenerator_.stateNames();
    const std::vector<std::string> &constantNames = cellmlSourceCodeGenerator_.constantNames();


    // -----------------------------
    // parse all entries of "mappings" to tupleEntries of the following form
    /*
    "mappings": {
      ("parameter",0): ("state","membrane/V"),
      ("parameter",1): ("algebraic","razumova/activestress"),
      ("parameter",2): ("constant",0),
      ("outputConnectorSlot",0): ("parameter",0),
      ("outputConnectorSlot",1): "razumova/l_hs",       #<-- this will be converted now
    }
    */
    using ValueTupleType = std::pair<std::string,PyObject *>;
    using KeyTupleType = std::pair<std::string,int>;

    std::vector<std::pair<KeyTupleType,ValueTupleType>> tupleEntries;

    std::pair<PyObject *,PyObject *> item = this->specificSettings_.template getOptionDictBegin<PyObject *,PyObject *>("mappings");

    for (;!this->specificSettings_.getOptionDictEnd("mappings");
         this->specificSettings_.template getOptionDictNext<PyObject *,PyObject *>("mappings",item))
    {
      KeyTupleType keyTuple;
      ValueTupleType valueTuple;

      if (PyTuple_Check(item.first))
      {
        keyTuple = PythonUtility::convertFromPython<KeyTupleType>::get(item.first);

        // if tuple is not parameter or outputConnectorSlot
        if (keyTuple.first != "parameter" && keyTuple.first != "outputConnectorSlot")
        {
          keyTuple.first = "parameter";
          LOG(ERROR) << this->specificSettings_ << "[\"mappings\"], key " << item.first << " is not a tuple (\"parameter\",index) "
            << "or (\"outputConnectorSlot\",index). Assuming (\"parameter\", " << keyTuple.second << ")";
        }
      }
      else
      {
        int index = PythonUtility::convertFromPython<int>::get(item.first);
        keyTuple = std::pair<std::string,int>("parameter", index);

        LOG(ERROR) << this->specificSettings_ << "[\"mappings\"], key " << item.first << " is not a tuple (\"parameter\",index) "
          << "or (\"outputConnectorSlot\",index). Assuming (\"parameter\", " << index << ")";
      }

      if (PyTuple_Check(item.second))
      {
        valueTuple = PythonUtility::convertFromPython<ValueTupleType>::get(item.second);
      }
      else
      {
        std::string name = PythonUtility::convertFromPython<std::string>::get(item.second);

        // find if the string is a algebraic name given in algebraicNames
        std::vector<std::string>::const_iterator pos = std::find(algebraicNames.begin(), algebraicNames.end(), name);
        if (pos != algebraicNames.end())
        {
          valueTuple.first = "algebraic";
          valueTuple.second = PythonUtility::convertToPython<std::string>::get(name);
        }
        else
        {
          std::vector<std::string>::const_iterator pos = std::find(stateNames.begin(), stateNames.end(), name);
          if (pos != stateNames.end())
          {
            valueTuple.first = "state";
            valueTuple.second = PythonUtility::convertToPython<std::string>::get(name);
          }
          else
          {
            std::vector<std::string>::const_iterator pos = std::find(constantNames.begin(), constantNames.end(), name);
            if (pos != constantNames.end())
            {
              valueTuple.first = "constant";
              valueTuple.second = PythonUtility::convertToPython<std::string>::get(name);
            }
            else
            {
              std::stringstream stateNamesListing;
              std::stringstream algebraicNamesListing;
              std::stringstream constantNamesListing;

              // compose message for state names
              const int nColumns = 3;
              int nPerColumn = std::ceil(double(stateNames.size())/nColumns);
              for (int i = 0; i < nPerColumn; i++)
              {
                for (int j = 0; j < nColumns; j++)
                {
                  int index = j*nPerColumn + i;
                  if (index < stateNames.size())
                    stateNamesListing << stateNames[index] << std::string(40-stateNames[index].length(), ' ');
                }
                stateNamesListing << "\n  ";
              }

              // compose message for algebraic names
              nPerColumn = std::ceil(double(algebraicNames.size())/nColumns);
              for (int i = 0; i < nPerColumn; i++)
              {
                for (int j = 0; j < nColumns; j++)
                {
                  int index = j*nPerColumn + i;
                  if (index < algebraicNames.size())
                    algebraicNamesListing << algebraicNames[index] << std::string(40-algebraicNames[index].length(), ' ');
                }
                algebraicNamesListing << "\n  ";
              }

              // compose message for constant names
              nPerColumn = std::ceil(double(constantNames.size())/nColumns);
              for (int i = 0; i < nPerColumn; i++)
              {
                for (int j = 0; j < nColumns; j++)
                {
                  int index = j*nPerColumn + i;
                  if (index < constantNames.size())
                    constantNamesListing << constantNames[index] << std::string(40-constantNames[index].length(), ' ');
                }
                constantNamesListing << "\n  ";
              }

              LOG(ERROR) << this->specificSettings_ << "[\"mappings\"], name \"" << name << "\" is neither state, algebraic nor constant.\n\n"
                << "  Valid state names are:\n  " << stateNamesListing.str() << "\n\n"
                << "  Valid algebraic names are:\n  " << algebraicNamesListing.str() << "\n\n"
                << "  Valid constant names are:\n  " << constantNamesListing.str() << "\n";
            }
          }
        }
      }

      tupleEntries.push_back(std::pair<KeyTupleType,ValueTupleType>(keyTuple, valueTuple));
    }

    // sort the tupleEntries by no of the key, i.e. parameter no. and output connector slot no.
    std::sort(tupleEntries.begin(), tupleEntries.end(), [](
      const std::pair<KeyTupleType,ValueTupleType> &a,
      const std::pair<KeyTupleType,ValueTupleType> &b)
    {
      return a.first.first.length() < b.first.first.length() ||
        (a.first.first.length() == b.first.first.length() && a.first.second < b.first.second);
    });

    LOG(DEBUG) << "sorted tupleEntries: " << tupleEntries;

    // -----------------------------
    // parse the strings in the entries and convert them to stateNo, algebraicNo or constantNo
    /*
    "mappings": {
      ("parameter",0): ("state","membrane/V"),                      #<-- this will be converted now
      ("parameter",1): ("algebraic","razumova/activestress"),    #<-- this will be converted now
      ("parameter",2): ("constant",0),
      ("outputConnectorSlot",0): ("parameter",0),
      ("outputConnectorSlot",1): ("algebraic", "razumova/l_hs"), #<-- this will be converted now
    }
    */
    // result is settings for parametersUsedAsAlgebraic, statesForTransfer, algebraicsForTransfer and parametersForTransfer
    parametersUsedAsAlgebraic.clear();
    parametersUsedAsConstant.clear();
    statesForTransfer.clear();
    algebraicsForTransfer.clear();
    parametersForTransfer.clear();

    for (std::pair<KeyTupleType,ValueTupleType> &tupleEntry : tupleEntries)
    {
      // parse value
      ValueTupleType valueTuple = tupleEntry.second;

      std::string specifier = valueTuple.first;    // state, algebraic, constant or parameter
      int fieldNo;              // no of the state, algebraic, constant or parameter

      // parse value tuple into (specifier and fieldNo)
      if (specifier == "state")
      {
        std::string stateName = PythonUtility::convertFromPython<std::string>::get(valueTuple.second);

        // if parsed string is a number, this is directly the state no
        if (!stateName.empty() && stateName.find_first_not_of("0123456789") == std::string::npos)
        {
          fieldNo = atoi(stateName.c_str());
        }
        else
        {
          // find if the string is a state name given in stateNames
          const std::vector<std::string>::const_iterator pos = std::find(stateNames.begin(), stateNames.end(), stateName);
          if (pos != stateNames.end())
          {
            fieldNo = std::distance(stateNames.begin(),pos);
          }
          else
          {
            LOG(ERROR) << "In " << this->specificSettings_ << "[\"mappings\"], state \"" << stateName << "\" is not valid. "
              << "Valid state names are: " << stateNames << ".\n"
              << "You can also specify the state no. (index of STATES[] array in C file) instead of the name.";
          }
        }
      }
      else if (specifier == "algebraic" || specifier == "intermediate")     // "intermediate" was the old name
      {
        std::string algebraicName = PythonUtility::convertFromPython<std::string>::get(valueTuple.second);

        // if parsed string is a number, this is directly the algebraic no
        if (!algebraicName.empty() && algebraicName.find_first_not_of("0123456789") == std::string::npos)
        {
          fieldNo = atoi(algebraicName.c_str());
        }
        else
        {
          // find if the string is a algebraic name given in algebraicNames
          std::vector<std::string>::const_iterator pos = std::find(algebraicNames.begin(), algebraicNames.end(), algebraicName);
          if (pos != algebraicNames.end())
          {
            fieldNo = std::distance(algebraicNames.begin(),pos);
          }
          else
          {
            LOG(ERROR) << "In " << this->specificSettings_ << "[\"mappings\"], algebraic \"" << algebraicName << "\" is not valid. "
              << "Valid algebraic names are: " << algebraicNames << ".\n"
              << "You can also specify the algebraic no. (index of ALGEBRAICS[] array in C file) instead of the name.";
          }
        }
      }
      else if (specifier == "constant")
      {
        std::string constantName = PythonUtility::convertFromPython<std::string>::get(valueTuple.second);

        // if parsed string is a number, this is directly the constant no
        if (!constantName.empty() && constantName.find_first_not_of("0123456789") == std::string::npos)
        {
          fieldNo = atoi(constantName.c_str());
        }
        else
        {
          // find if the string is a constant name given in constantNames
          std::vector<std::string>::const_iterator pos = std::find(constantNames.begin(), constantNames.end(), constantName);
          if (pos != constantNames.end())
          {
            fieldNo = std::distance(constantNames.begin(),pos);
          }
          else
          {
            LOG(ERROR) << "In " << this->specificSettings_ << "[\"mappings\"], constant \"" << constantName << "\" is not valid. "
              << "Valid constant names are: " << constantNames << ".\n"
              << "You can also specify the constant no. (index of CONSTANTS[] array in C file) instead of the name.";
          }
        }
      }
      else
      {
        LOG(FATAL) << "In " << this->specificSettings_ << "[\"mappings\"], invalid field \"" << specifier << "\", allowed values are: "
          << "\"state\", \"algebraic\", \"constant\".";
      }

      // store according to key
      KeyTupleType keyTuple = tupleEntry.first;

      // handle parameter
      if (keyTuple.first == "parameter")
      {
        int parameterNo = keyTuple.second;    // the no. of the parameter

        LOG(DEBUG) << "parameter " << parameterNo << " to \"" << valueTuple.first << "\", fieldNo " << fieldNo;

        if (valueTuple.first == "algebraic" || valueTuple.first == "intermediate")
        {
          parametersUsedAsAlgebraic.push_back(fieldNo);
        }
        else if (valueTuple.first == "constant")
        {
          parametersUsedAsConstant.push_back(fieldNo);
        }
        else
        {
          LOG(ERROR) << "In " << this->specificSettings_ << "[\"mappings\"], you can map "
            << "paramaters only to algebraics and constants, not states.";
        }
      }
      else
      {
        // output connector slots

        int slotNo = keyTuple.second;    // the no. of the output connector slot

        LOG(DEBUG) << "output connector slot " << slotNo << " to \"" << valueTuple.first << "\", fieldNo " << fieldNo;

        if (valueTuple.first == "state")
        {
          statesForTransfer.push_back(fieldNo);
        }
        else if (valueTuple.first == "algebraic" || valueTuple.first == "intermediate")
        {
          algebraicsForTransfer.push_back(fieldNo);
        }
        else if (valueTuple.first == "parameter")
        {
          parametersForTransfer.push_back(fieldNo);
        }
        else
        {
          std::vector<int>::iterator parameter = std::find(parametersUsedAsConstant.begin(), parametersUsedAsConstant.end(), fieldNo);

          if (parameter != parametersUsedAsConstant.end())
          {
            int parameterNo = parametersUsedAsAlgebraic.size() + parameter - parametersUsedAsConstant.begin();
            parametersForTransfer.push_back(parameterNo);
            LOG(DEBUG) << "Constant " << valueTuple.second << " was set as output connector slot " << slotNo
              << " and is mapped to parameter " << parameterNo;
          }
          else
          {
            LOG(ERROR) << "In " << this->specificSettings_ << "[\"mappings\"], you can only connect states, algebraics and parameters "
              << "to outputConnectorSlots, not constants (\"" << valueTuple.second << "\" is a " << valueTuple.first << "). \n"
              << "You can use it for a slot if you define it as parameter first.";
          }
        }
      }
    }

    // create a string of parameter and contstant mappings
    std::stringstream s;
    for (int i = 0; i < parametersUsedAsAlgebraic.size(); i++)
    {
      int algebraicNo = parametersUsedAsAlgebraic[i];
      if (algebraicNo != -1 && algebraicNo < nAlgebraics_)
      {
        s << "  parameter " << i << " maps to algebraic " << algebraicNo << " (\"" << algebraicNames[algebraicNo] << "\")\n";
        parameterNames.push_back(algebraicNames[algebraicNo]);
      }
    }
    for (int i = 0; i < parametersUsedAsConstant.size(); i++)
    {
      int constantNo = parametersUsedAsConstant[i];
      if (constantNo != -1 && constantNo < constantNames.size())
      {
        s << "  parameter " << parametersUsedAsAlgebraic.size()+i << " maps to constant " << constantNo << " (\"" << constantNames[constantNo] << "\")\n";
        parameterNames.push_back(constantNames[constantNo]);
      }
    }

    LOG(DEBUG) << s.str();
  }
  else
  {
    this->specificSettings_.getOptionVector("parametersUsedAsAlgebraic", parametersUsedAsAlgebraic);
    this->specificSettings_.getOptionVector("parametersUsedAsConstant", parametersUsedAsConstant);

    this->specificSettings_.template getOptionVector<int>("statesForTransfer", statesForTransfer);
    this->specificSettings_.template getOptionVector<int>("algebraicsForTransfer", algebraicsForTransfer);
    this->specificSettings_.template getOptionVector<int>("parametersForTransfer", parametersForTransfer);
  }

  LOG(DEBUG) << "parametersUsedAsAlgebraic: " << parametersUsedAsAlgebraic;
  LOG(DEBUG) << "parametersUsedAsConstant:  " << parametersUsedAsConstant;
  LOG(DEBUG) << "statesForTransfer:         " << statesForTransfer;
  LOG(DEBUG) << "algebraicsForTransfer:     " << algebraicsForTransfer;
  LOG(DEBUG) << "parametersForTransfer:     " << parametersForTransfer;

  if (this->specificSettings_.hasKey("outputAlgebraicIndex") || this->specificSettings_.hasKey("outputIntermediateIndex"))
  {
    LOG(WARNING) << specificSettings_ << "[\"outputAlgebraicIndex\"] is no longer a valid option, use \"algebraicsForTransfer\" instead!";
  }

  if (this->specificSettings_.hasKey("outputStateIndex"))
  {
    LOG(WARNING) << specificSettings_ << "[\"outputStateIndex\"] is no longer a valid option, use \"statesForTransfer\" instead!";
  }
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
template<typename FunctionSpaceType2>
bool CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nStates_>> initialValues)
{
  LOG(TRACE) << "CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::setInitialValues";

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

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
functionSpace()
{
  return functionSpace_;
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
getNumbers(int& nInstances, int& nAlgebraics, int& nParameters)
{
  nInstances = nInstances_;
  nAlgebraics = nAlgebraics_;
  nParameters = cellmlSourceCodeGenerator_.nParameters();
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
std::vector<int> &CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
statesForTransfer()
{
  return data_.statesForTransfer();
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
std::vector<int> &CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
algebraicsForTransfer()
{
  return data_.algebraicsForTransfer();
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
getStateNames(std::vector<std::string> &stateNames)
{
  stateNames = this->cellmlSourceCodeGenerator_.stateNames();
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
typename CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::Data &CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
data()
{
  return this->data_;
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
CellmlSourceCodeGenerator &CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
cellmlSourceCodeGenerator()
{
  return this->cellmlSourceCodeGenerator_;
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
PythonConfig CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
specificSettings()
{
  return this->specificSettings_;
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
constexpr int  CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
nAlgebraics() const
{
  return nAlgebraics_;
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
std::shared_ptr<typename CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::OutputConnectorDataType>
CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::
getOutputConnectorData()
{
  return this->data_.getOutputConnectorData();
}
