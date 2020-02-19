#pragma once

#include <Python.h>  // has to be the first included header

#include <vector>

#include "control/dihu_context.h"
#include "output_writer/manager.h"
#include "function_space/function_space.h"
#include "data_management/cellml_adapter.h"
#include "data_management/output_connector_data.h"
#include "cellml/source_code_generator/source_code_generator.h"


/** This is the base class of the CellmlAdapter, that handles common functionality.
 * nStates_: number of states in one instance of the CellML problem
 * 
 *  Naming:
 *   Intermediate (opendihu) = KNOWN (OpenCMISS) = Algebraic (OpenCOR)
 *   Parameter (opendihu, OpenCMISS) = KNOWN (OpenCMISS), in OpenCOR also algebraic
 *   Constant - these are constants that are only present in the source files
 *   State: state variable
 *   Rate: the time derivative of the state variable, i.e. the increment value in an explicit Euler stepping
 */
template <int nStates_, int nIntermediates_, typename FunctionSpaceType>
class CellmlAdapterBase
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,nIntermediates_> FieldVariableIntermediates;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nStates_> FieldVariableStates;
  typedef Data::CellmlAdapter<nStates_, nIntermediates_, FunctionSpaceType> Data;

/** The data type of the output connector of the CellML adapter.
 *  This is the data that will be transferred to connected solvers.
 *  The first value, value0, is the state variable and can, e.g., be configured to contain Vm (by setting "outputStateIndex" in python settings).
 *  The second value, value1, can, e.g., be configured to contain alpha (by setting "outputIntermediateIndex" in python settings).
 */
  typedef ::Data::OutputConnectorData<FunctionSpaceType,nStates_,nIntermediates_> OutputConnectorDataType;

  //! constructor from context
  CellmlAdapterBase(DihuContext context, bool noNewOutputWriter);

  //! constructor from context
  CellmlAdapterBase(DihuContext context);

  //! destructor
  ~CellmlAdapterBase();

  //! return the compile-time constant number of state variables of one instance that will be integrated
  static constexpr int nComponents();

  //! return the compile-time constant number of state variables of one instance that will be integrated
  static constexpr int nStates();

  //! load model, use settings given in context
  void initialize();

  //! create generic function space with given number of instances
  void initializeFromNInstances(int nInstances);
  
  //! set initial values as given in python config
  template<typename FunctionSpaceType2>
  bool setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nStates_>> initialValues);

  //! set the solution field variable in the data object, that actual data is stored in the timestepping scheme object
  void setSolutionVariable(std::shared_ptr<FieldVariableStates> states);

  //! pass on the output connector data object from the timestepping scheme object to be modified,
  //! if there are intermediates for transfer, they will be set in the outputConnectorDataTimeStepping
  void setOutputConnectorData(std::shared_ptr<::Data::OutputConnectorData<FunctionSpaceType,nStates_>> outputConnectorDataTimeStepping);

  //! return the mesh
  std::shared_ptr<FunctionSpaceType> functionSpace();

  //! get number of instances, number of intermediates and number of parameters
  void getNumbers(int &nInstances, int &nIntermediates, int &nParameters);

  //! return references to statesForTransfer and intermediatesForTransfer, the states and intermediates that should be used for output connector data transfer
  void getStatesIntermediatesForTransfer(std::vector<int> &statesForTransfer, std::vector<int> &intermediatesForTransfer);

  //! get a vector with the names of the states
  void getStateNames(std::vector<std::string> &stateNames);

  //! get the const number of intermediates
  constexpr int nIntermediates() const;

  //! return reference to the data object that stores all field variables
  Data &data();

  //! get the python config object that contains all python settings for the CellML adapter
  PythonConfig specificSettings();

  //! get the initialized cellmlSourceCodeGenerator
  CellmlSourceCodeGenerator &cellmlSourceCodeGenerator();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the output_connector_data class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

protected:

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer

  std::shared_ptr<FunctionSpaceType> functionSpace_;    ///< a mesh, there are as many instances of the same CellML problem as there are nodes in the mesh
  Data data_;     ///< the data object that stores all variables, i.e. intermediates and states

  int nInstances_;         ///< number of instances of the CellML problem. Usually it is the number of mesh nodes when a mesh is used. When running in parallel this is the local number of instances without ghosts.

  int internalTimeStepNo_ = 0; ///< the counter how often the right hand side was called

  CellmlSourceCodeGenerator cellmlSourceCodeGenerator_;    ///< object that holds all source code related to the model
};

#include "cellml/00_cellml_adapter_base.tpp"
