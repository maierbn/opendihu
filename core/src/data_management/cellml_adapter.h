#pragma once

#include <Python.h>  // has to be the first included header

#include "slot_connection/slot_connector_data.h"
#include "control/python_config/python_config.h"
#include "data_management/data.h"

namespace Data
{

/**  The datastructures used for the cellml adapter.
 *   A DiscretizableInTime class has no own slotConnectorData, since is has no own data.
 *   The relevant data (solution vector) is stored in the timestepping scheme. However, the algebraics and parameters are stored here (in the data object).
 *   Therefore, the CellmlAdapter has a copy of the slotConnectorData of the timestepping scheme
 */
template <int nStates, int nAlgebraics, typename FunctionSpaceType>
class CellmlAdapter : public Data<FunctionSpaceType>
{
public:
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nAlgebraics> FieldVariableAlgebraics;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nStates> FieldVariableStates;
  //typedef SlotConnectorData<FunctionSpaceType,nStates,nAlgebraics> SlotConnectorDataType;

  //! constructor
  CellmlAdapter(DihuContext context);

  //! constructor as a copy
  CellmlAdapter(DihuContext context, const CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType> &rhs);

  //! return a reference to the algebraics field variable
  std::shared_ptr<FieldVariableAlgebraics> algebraics();

  //! return the field variable that stores all states
  std::shared_ptr<FieldVariableStates> states();

  //! return the field variable that stores all parameters
  std::shared_ptr<FieldVariableAlgebraics> parameters();

  //! get a raw r/w memory pointer to the values of the parameters field variable, to be used in the rhs function
  //! the array pointed to by this pointer contains parameter values for all instances, nAlgebraics parameters for each instance, in struct-of-array memory layout
  double *parameterValues();

  //! initialize data structures
  void initialize();

  //! assign the field variable that stores the states, this variable is owned by the timestepping scheme
  void setStatesVariable(std::shared_ptr<FieldVariableStates> states);

  //! give the names of all algebraics, will be called before initialize()
  void setAlgebraicAndParameterNames(const std::vector<std::string> &algebraicNames, const std::vector<std::string> &parameterNames, const std::vector<std::string> &slotNames);

  //! get the parameteValues_ pointer from the parameters field variable, then the field variable can no longer be used until restoreParameterValues() gets called
  void prepareParameterValues();

  //! restore the parameterValues_ pointer, such that the field variable can be used again
  void restoreParameterValues();

  //! return a reference to statesForTransfer_
  std::vector<int> &statesForTransfer();

  //! return a reference algebraicsForTransfer_
  std::vector<int> &algebraicsForTransfer();

  //! return a reference parametersForTransfer_
  std::vector<int> &parametersForTransfer();

  //! initialize the slotConnectorData object of the timestepping scheme, this method is called once in the initialize() of the timestepping scheme
  //! It adds all slots that were defined in the CellmlAdapter to the slotConnectorData of the timestepping scheme
  void setSlotConnectorDataTimestepping(std::shared_ptr<SlotConnectorData<FunctionSpaceType,nStates>> slotConnectorDataTimeStepping);

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,     // geometry
    std::shared_ptr<FieldVariableAlgebraics>,     // algebraics
    std::shared_ptr<FieldVariableStates>,            // states
    std::shared_ptr<FieldVariableAlgebraics>          // parameters
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  //! initialize the slotConnectorData_ object, this can only be done after the states variable has been set
  void initializeSlotConnectorData();
  
  //! retrieve the current parameters_ field variable from the slot connector, because is could have been changed by the slot connector transfer
  void updateParametersPointer();

  std::shared_ptr<FieldVariableAlgebraics> algebraics_;   //< algebraics field variable
  std::shared_ptr<FieldVariableStates> states_;           //< states field variable, this is a shared pointer with the timestepping scheme, which own the actual variable (creates it)
  std::shared_ptr<FieldVariableAlgebraics> parameters_;   //< parameters field variable, the number of components is equal or less than the number of algebraics in order to not have to specify the number of parameters at compile time. This possibly creates a vector that is too large which is not harmful.
  double *parameterValues_;                               //< a pointer to the data of the parameters_ Petsc Vec of the field variable
  std::vector<std::string> algebraicNames_;               //< component names of the algebraics field variable
  std::vector<std::string> parameterNames_;               //< component names of the parameter field variable
  std::vector<std::string> slotNames_;                    //< names of the data slots that are used for slot connectors

  std::shared_ptr<SlotConnectorData<FunctionSpaceType,nStates>> slotConnectorDataTimestepping_;   //< slotConnectorData = the object that holds all components of field variables that will be transferred to other solvers, this is a copy of the slot connector data object of the timestepping scheme
  PythonConfig specificSettings_;                               //< the settings object

  std::vector<int> statesForTransfer_;                    //< state no.s to transfer to other solvers within slot connector data
  std::vector<int> algebraicsForTransfer_;                //< algebraic no.s to transfer to other solvers within slot connector data
  std::vector<int> parametersForTransfer_;                //< parameter no.s to transfer to other solvers within slot connector data
};

} // namespace Data

#include "data_management/cellml_adapter.tpp"
