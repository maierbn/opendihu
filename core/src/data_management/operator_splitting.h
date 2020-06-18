#pragma once

#include <Python.h>  // has to be the first included header
#include <tuple>

#include "control/dihu_context.h"

namespace Data
{

/** The data class used by operator splitting and coupling.
 *  This class only manages the outputConnectorData, besides that
 *  the class does not contain anything.
  */
template<typename TimeStepping1, typename TimeStepping2>
class OperatorSplitting :
  public Data<typename TimeStepping1::FunctionSpace>
{
public:
  typedef std::tuple<
    std::shared_ptr<typename TimeStepping1::OutputConnectorDataType>,
    std::shared_ptr<typename TimeStepping2::OutputConnectorDataType>
  > OutputConnectorDataType;

  //! constructor
  OperatorSplitting(DihuContext context);

  //! initiialize the output connector data
  void initialize(TimeStepping1 &timeStepping1, TimeStepping2 &timeStepping2);

  //! this method is needed but does nothing here, since there are no field variables in this data object
  void createPetscObjects(){}

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

protected:

  std::shared_ptr<OutputConnectorDataType> outputConnectorData_;  //< the object that holds output connector data that will be transferred between solvers
};

} // namespace Data

#include "data_management/operator_splitting.tpp"
