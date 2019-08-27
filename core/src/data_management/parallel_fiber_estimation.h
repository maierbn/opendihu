#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "field_variable/field_variable.h"
#include "quadrature/gauss.h"

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

  typedef SpatialDiscretization::FiniteElementMethod<
    typename FunctionSpaceType::Mesh,
    typename FunctionSpaceType::BasisFunction,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > FiniteElementMethodType;

  //! constructor
  ParallelFiberEstimation(DihuContext context);

  //! destructur
  ~ParallelFiberEstimation();

  //! return a reference to the gradient vector field
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> gradient();

  //! return a reference to the dirichletValues field
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> dirichletValues();

  //! return a reference to the field for the condition number of the jacobian
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> jacobianConditionNumber();

  //! set the problem variable
  void setProblem(std::shared_ptr<FiniteElementMethodType> problem);

  //! print all stored data to stdout
  virtual void print();

  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> FunctionSpaceFiber;
  typedef FieldVariable::FieldVariable<FunctionSpaceFiber,3> FieldVariableFiberGeometry;
  
  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,  // geometry
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>,  // solution
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>,   // rhs
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>,   // rhs neumann bc
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,  // gradient field
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>,  // dirichlet values
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>  // dirichlet values
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

protected:

  //! initializes the vectors with size
  virtual void createPetscObjects();

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> gradient_;  ///< the gradient field of the Laplace flow solution
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> dirichletValues_;  ///< values of dirichlet BC or -1, where no dirichlet BC is prescribed
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> jacobianConditionNumber_;  ///< condition number of the jacobian at each point

  std::shared_ptr<FiniteElementMethodType> problem_;   ///< the DiscretizableInTime object that is used for FE solution

};

} // namespace Data

#include "data_management/parallel_fiber_estimation.tpp"
