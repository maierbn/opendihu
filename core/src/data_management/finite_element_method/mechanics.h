#pragma once

#include <Python.h>  // has to be the first included header

namespace Data
{

/**  The datastructures used for mechanics solver
  */
template<typename FunctionSpaceType>
class Mechanics : public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,FunctionSpaceType::dim()> FieldVariableType;

  //! constructor
  Mechanics(DihuContext context);

  //! return a reference to the displacements vector
  std::shared_ptr<FieldVariableType> displacements();

  //! return a reference to the rhs
  std::shared_ptr<FieldVariableType> rightHandSide();

  //! get a reference to the rhs matrix
  Mat &stiffnessMatrix();

  //! initialize
  void initialize();

  //! reset to uninitialized state
  void reset();

  //! print all stored data to stdout
  void print();

  //! set the external forces rhs ∫ (T_L*phi_L)_a δu_aM phi_M ds,   a = 1,2,3, L,M ∈ surface dofs, that was computed by NeumannBoundaryConditions
  //! this sets the negative value in rightHandSide_
  void setExternalForcesRightHandSide(std::shared_ptr<FieldVariableType> rhs);

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,  // geometry
    std::shared_ptr<FieldVariableType>,  // displacements
    std::shared_ptr<FieldVariableType>   // rhs
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  Mat stiffnessMatrix_;   ///< the stiffness matrix, compound of submatrices for the D dimensions (MatNest)
  std::array<std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>>, FunctionSpaceType::dim()*FunctionSpaceType::dim()> stiffnessMatrixComponents_; ///< submatrices in stiffnessMatrix

  std::shared_ptr<FieldVariableType> displacements_; ///< solution of the mechanics problem
  std::shared_ptr<FieldVariableType> rightHandSide_; ///< the right hand side which contains the external forces
};

} // namespace Data

#include "data_management/mechanics.tpp"
