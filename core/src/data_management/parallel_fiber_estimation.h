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

/**  The datastructures used for streamline tracer.
 *   BaseDataType is a Data class that provides the solution field variable for the streamline tracer to operate on.
 */
template<typename FunctionSpaceType>
class ParallelFiberEstimation :
  public Data<FunctionSpaceType>
{
public:

  //! constructor
  ParallelFiberEstimation(DihuContext context);

  //! destructur
  ~ParallelFiberEstimation();

  //! return a reference to the solution vector, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> gradient();

  //! print all stored data to stdout
  virtual void print();

  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> FunctionSpaceFibre;
  typedef FieldVariable::FieldVariable<FunctionSpaceFibre,3> FieldVariableFibreGeometry;
  
  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>  // gradient field
  > OutputFieldVariables;


  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

protected:

  //! initializes the vectors with size
  virtual void createPetscObjects();

};

} // namespace Data

#include "data_management/parallel_fiber_estimation.tpp"
