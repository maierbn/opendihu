#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "output_connector_data_transfer/output_connector_data_transfer.h"
#include "data_management/finite_element_method/finite_elements.h"

/** Transfer between one scalar field variable (e.g. from finite element method) and one field variables with given component number and prefactor (e.g. from cellml)
 */
template<typename FunctionSpaceType1, typename OutputConnectorDataType2>
class SolutionVectorMapping<
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>,   // <fieldVariableType,componentNo>
  OutputConnectorDataType2
> :
  public SolutionVectorMapping<
    Data::FieldVariableComponent<FunctionSpaceType1,1>,
    OutputConnectorDataType2
  >
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> &transferableSolutionData1,
                       const OutputConnectorDataType2 &transferableSolutionData2,
                       const std::string transferSlotName);
};

/** Transfer between one scalar field variable (e.g. from finite element method) and one field variables with given component number and prefactor (e.g. from cellml)
 *  The other way
 */
template<typename OutputConnectorDataType1, typename FunctionSpaceType2>
class SolutionVectorMapping<
  OutputConnectorDataType1,
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>    // <fieldVariableType,componentNo>
> :
  public SolutionVectorMapping<
    OutputConnectorDataType1,
    Data::FieldVariableComponent<FunctionSpaceType2,1>
  >
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const OutputConnectorDataType1 &transferableSolutionData1,
                       const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>> &transferableSolutionData2,
                       const std::string transferSlotName);

};

/** Transfer between one scalar field variable (e.g. from finite element method) and one field variables with given component number and prefactor (e.g. from cellml)
 *  The other way
 */
template<typename FunctionSpaceType1, typename FunctionSpaceType2>
class SolutionVectorMapping<
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>,
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>
> :
  public SolutionVectorMapping<
    Data::FieldVariableComponent<FunctionSpaceType1,1>,
    Data::FieldVariableComponent<FunctionSpaceType2,1>
  >
{
public:
  //! transfer the data from transferableSolutionData1 to transferableSolutionData2, as efficient as possible, where there are multiple slots that could be transferred (e.g. at cellmlAdapter), use the one specified by transferSlotName
  static void transfer(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> &transferableSolutionData1,
                       const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>> &transferableSolutionData2,
                       const std::string transferSlotName);

};

#include "output_connector_data_transfer/output_connector_data_transfer_finite_element_method_cellml.tpp"
