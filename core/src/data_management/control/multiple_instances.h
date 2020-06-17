#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  Multiple instances of a time stepping
 */
template<typename FunctionSpaceType, typename BaseTimesteppingType>
class MultipleInstances :
  public Data<FunctionSpaceType>
{
public:
  typedef typename BaseTimesteppingType::Data BaseDataType;

  //! constructor
  MultipleInstances(DihuContext context);

  //! destructur
  ~MultipleInstances();

  //! print all stored data to stdout
  virtual void print();

  //! set the data objects of the instances
  void setInstancesData(std::vector<BaseTimesteppingType> &instances);

  //! field variables that will be output by outputWriters
  typedef std::tuple<std::vector<typename BaseDataType::FieldVariablesForOutputWriter>> FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

protected:

  //! initializes the vectors with size
  virtual void createPetscObjects();

  std::vector<std::shared_ptr<BaseDataType>> instancesData_;    //< the data objects of all instances
};

} // namespace Data

#include "data_management/control/multiple_instances.tpp"
