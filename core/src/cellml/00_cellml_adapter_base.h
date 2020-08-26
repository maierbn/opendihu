#pragma once

#include <Python.h>  // has to be the first included header

#include <vector>

#include "control/dihu_context.h"
#include "interfaces/discretizable_in_time.h"
#include "output_writer/manager.h"
#include "function_space/function_space.h"
#include "data_management/cellml_adapter.h"
#include "slot_connection/slot_connector_data.h"
#include "cellml/source_code_generator/source_code_generator.h"


/** This is the base class of the CellmlAdapter, that handles common functionality.
 * nStates_: number of states in one instance of the CellML problem
 * 
 *  Naming:
 *   Algebraic (opendihu) = KNOWN (OpenCMISS) = Algebraic (OpenCOR)
 *   Parameter (opendihu, OpenCMISS) = KNOWN (OpenCMISS), in OpenCOR also algebraic
 *   Constant - these are constants that are only present in the source files
 *   State: state variable
 *   Rate: the time derivative of the state variable, i.e. the increment value in an explicit Euler stepping
 */
template <int nStates_, int nAlgebraics_, typename FunctionSpaceType>
class CellmlAdapterBase :
  public DiscretizableInTime<FunctionSpaceType,nStates_>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,nAlgebraics_> FieldVariableAlgebraics;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nStates_> FieldVariableStates;
  typedef Data::CellmlAdapter<nStates_, nAlgebraics_, FunctionSpaceType> Data;

/** The data type of the slot connector of the CellML adapter.
 *  This is the data that will be transferred to connected solvers.
 *  The first value, value0, is the state variable and can, e.g., be configured to contain Vm (by setting "outputStateIndex" in python settings).
 *  The second value, value1, can, e.g., be configured to contain alpha (by setting "outputAlgebraicIndex" in python settings).
 */
  typedef ::Data::SlotConnectorData<FunctionSpaceType,nStates_,nAlgebraics_> SlotConnectorDataType;

  //! constructor from context
  CellmlAdapterBase(DihuContext context);

  //! constructor with given data
  CellmlAdapterBase(DihuContext context, const CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::Data &rhsData);

  //! destructor
  ~CellmlAdapterBase();

  //! return the compile-time constant number of state variables of one instance that will be integrated
  static constexpr int nComponents();

  //! return the compile-time constant number of state variables of one instance that will be integrated
  static constexpr int nStates();

  //! load model, use settings given in context
  void initialize();

  //! set initial values as given in python config
  bool setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStates_>> initialValues);

  //! initialize all information from python settings key "mappings", this sets parametersUsedAsAlgebraics/States and outputAlgebraic/StatesIndex
  void initializeMappings(std::vector<int> &parametersUsedAsAlgebraic, std::vector<int> &parametersUsedAsConstant,
                          std::vector<int> &statesForTransfer, std::vector<int> &algebraicsForTransfer, std::vector<int> &parametersForTransfer,
                          std::vector<std::string> &parameterNames, std::vector<std::string> &slotNames);

  //! set the solution field variable in the data object, that actual data is stored in the timestepping scheme object
  void setSolutionVariable(std::shared_ptr<FieldVariableStates> states);

  //! pass on the slot connector data object from the timestepping scheme object to be modified,
  //! if there are algebraics for transfer, they will be set in the slotConnectorDataTimeStepping
  void setSlotConnectorData(std::shared_ptr<::Data::SlotConnectorData<FunctionSpaceType,nStates_>> slotConnectorDataTimeStepping);

  //! return the mesh
  std::shared_ptr<FunctionSpaceType> functionSpace();

  //! get number of instances, number of algebraics and number of parameters
  void getNumbers(int &nInstances, int &nAlgebraics, int &nParameters);

  //! return a reference to statesForTransfer, the states that should be used for slot connector data transfer
  std::vector<int> &statesForTransfer();

  //! return a reference to algebraicsForTransfer, the algebraics that should be used for slot connector data transfer
  std::vector<int> &algebraicsForTransfer();

  //! get a vector with the names of the states
  void getStateNames(std::vector<std::string> &stateNames);

  //! get the const number of algebraics
  constexpr int nAlgebraics() const;

  //! return reference to the data object that stores all field variables
  Data &data();

  //! get the python config object that contains all python settings for the CellML adapter
  PythonConfig specificSettings();

  //! get the initialized cellmlSourceCodeGenerator
  CellmlSourceCodeGenerator &cellmlSourceCodeGenerator();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the slot_connector_data class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

protected:

  //! compute equilibrium of states for option "initializeStatesToEquilibrium"
  virtual void initializeToEquilibriumValues(std::array<double,nStates_> &statesInitialValues) = 0;


  DihuContext context_;                                    //< object that contains the python config for the current context and the global singletons meshManager and solverManager
  PythonConfig specificSettings_;                          //< python object containing the value of the python config dict with corresponding key
  OutputWriter::Manager outputWriterManager_;              //< manager object holding all output writer

  std::shared_ptr<FunctionSpaceType> functionSpace_;       //< a mesh, there are as many instances of the same CellML problem as there are nodes in the mesh
  Data data_;                                              //< the data object that stores all variables, i.e. algebraics and states
  static std::array<double,nStates_> statesInitialValues_; //< the initial values for the states, see setInitialValues
  static bool statesInitialValuesInitialized_;             //< if the statesInitialValues_ variables has been initialized

  int nInstances_;                                         //< number of instances of the CellML problem. Usually it is the number of mesh nodes when a mesh is used. When running in parallel this is the local number of instances without ghosts.
  int internalTimeStepNo_ = 0;                             //< the counter how often the right hand side was called

  bool initializeStatesToEquilibrium_;                     //< if the initial states should be computed until they reach an equilibrium
  double initializeStatesToEquilibriumTimestepWidth_;      //< timestep width for computation of equilibrium states

  CellmlSourceCodeGenerator cellmlSourceCodeGenerator_;    //< object that holds all source code related to the model
  bool initialized_;                                       //< if initialize has been called
};

#include "cellml/00_cellml_adapter_base.tpp"
