#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "field_variable/field_variable.h"
#include "function_space/function_space.h"
#include "data_management/output_connector_data.h"

namespace Data
{

/**  The datastructures for the febio adapter where febio computes a static hyperelastic problem.
 */
class QuasiStaticNonlinearElasticityFebio :
  public Data<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>>>
{
public:
  typedef ::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> FunctionSpace;

  typedef FieldVariable::FieldVariable<FunctionSpace,1> FieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpace,3> FieldVariableTypeVector;
  typedef FieldVariable::FieldVariable<FunctionSpace,6> FieldVariableTypeTensor;

  typedef OutputConnectorData<FunctionSpace,1> OutputConnectorDataType;

  //! constructor
  QuasiStaticNonlinearElasticityFebio(DihuContext context);

  //! return the field variable of the activation factor
  std::shared_ptr<FieldVariableType> activation();

  //! return the field variable
  std::shared_ptr<FieldVariableTypeVector> referenceGeometry();

  //! return the field variable
  std::shared_ptr<FieldVariableTypeVector> displacements();

  //! return the field variable
  std::shared_ptr<FieldVariableTypeVector> reactionForce();

  //! return the field variable
  std::shared_ptr<FieldVariableTypeTensor> cauchyStress();

  //! return the field variable
  std::shared_ptr<FieldVariableTypeTensor> pk2Stress();

  //! return the field variable
  std::shared_ptr<FieldVariableTypeTensor> greenLagrangeStrain();

  //! return the field variable
  std::shared_ptr<FieldVariableType> relativeVolume();

  //! initialize
  void initialize();

  //! print all stored data to stdout
  void print();

  //! get the output connection da
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
      std::shared_ptr<FieldVariableTypeVector>,        // geometry field
      std::shared_ptr<FieldVariableType>,              // activation
      std::shared_ptr<FieldVariableTypeVector>,        // displacements
      std::shared_ptr<FieldVariableTypeVector>,        // reactionForce
      std::shared_ptr<FieldVariableTypeTensor>,        // cauchyStress
      std::shared_ptr<FieldVariableTypeTensor>,        // PK2-Stress
      std::shared_ptr<FieldVariableTypeTensor>,        // greenLagrangeStrain
      std::shared_ptr<FieldVariableType>               // relativeVolume
    >
   FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<FieldVariableType> activation_;                 //< field variable of the activation factor field
  std::shared_ptr<FieldVariableTypeVector> referenceGeometry_;    //< field variable of the geometry value
  std::shared_ptr<FieldVariableTypeVector> displacements_;        //< field variable of the displacements
  std::shared_ptr<FieldVariableTypeVector> reactionForce_;        //< field variable of the reaction forces
  std::shared_ptr<FieldVariableTypeTensor> pk2Stress_;            //< field variable of the 2nd Piola-Kirchhoff stress, S
  std::shared_ptr<FieldVariableTypeTensor> cauchyStress_;         //< field variable of the Cauchy stress, sigma
  std::shared_ptr<FieldVariableTypeTensor> greenLagrangeStrain_;  //< field variable of the Green-Lagrange strain, E
  std::shared_ptr<FieldVariableType> relativeVolume_;             //< field variable of the relative volume (determinant of deformation gradient)

  std::shared_ptr<OutputConnectorDataType> outputConnectorData_;  //< the object that stores all components of field variables that will be transferred to other solvers
};

} // namespace Data

