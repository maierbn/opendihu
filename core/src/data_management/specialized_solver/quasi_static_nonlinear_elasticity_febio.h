#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "field_variable/field_variable.h"
#include "function_space/function_space.h"

namespace Data
{

/**  The datastructures for the febio adapter where febio computes a static hyperelastic problem.
 */
class QuasiStaticNonlinearElasticityFebio :
  public Data<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>>>
{
public:
  typedef ::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>> FunctionSpace;

  typedef FieldVariable::FieldVariable<FunctionSpace,1> FieldVariableType;

  //! constructor
  QuasiStaticNonlinearElasticityFebio(DihuContext context);

  //! return the field variable of the activation factor
  std::shared_ptr<FieldVariableType> activation();

  //! initialize
  void initialize();

  //! print all stored data to stdout
  void print();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
      std::shared_ptr<FieldVariableType>              // activation
    >
   FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<FieldVariableType> activation_; ///< field variable of the activation factor field

};

} // namespace Data

