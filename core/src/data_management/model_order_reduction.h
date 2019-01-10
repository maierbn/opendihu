#pragma once

#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"
#include "field_variable/field_variable.h"
#include "function_space/function_space.h"


namespace Data{

template<typename FunctionSpaceRowsType>  
class ModelOrderReduction:
  public Data<FunctionSpace::Generic>
{
public:
  typedef FieldVariable::FieldVariable<::FunctionSpace::Generic,1> FieldVariableType;
  
  //! constructor
  ModelOrderReduction(DihuContext context);
   
  virtual ~ModelOrderReduction();
  
  //! initialize the function space for rows of the snapshot matrix
  virtual void setFunctionSpaceRows(std::shared_ptr<FunctionSpaceRowsType> functionSpace);
   
  //! Basis for the reduced solution, V
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceRowsType,::FunctionSpace::Generic>> &basis();
   
  //! Transpose of the basis, V ^T
  std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,FunctionSpaceRowsType>> &basisTransp(); 
   
  //! initializes the basis V from an already existant Petsc Mat !?
  //void initializeBasis(Mat &basis);
   
  //! Reduced system matrix, A_R
  std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,::FunctionSpace::Generic>> &redSystemMatrix();
   
  //! initializes the redSysMatrix from an already existant Petsc Mat !?
  //void initializeRedSysMatrix(Mat &A_R);
   
  //! the reduced solution
  std::shared_ptr<FieldVariableType> &redSolution();
   
  //! The reduced order increment
  std::shared_ptr<FieldVariableType> &redIncrement();
   
  virtual void initialize() override;
   
private:
   
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceRowsType,::FunctionSpace::Generic>> basis_; // V
  std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,FunctionSpaceRowsType>> basisTransp_; // V^T
  std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,::FunctionSpace::Generic>> redSystemMatrix_;
   
  std::shared_ptr<FieldVariableType> redSolution_; //reduced solution
  std::shared_ptr<FieldVariableType> redIncrement_; //reduced increment
  
  std::shared_ptr<FunctionSpaceRowsType> functionSpaceRows_;
  
  //! Create the matrices and vectors for model order reduction
  void createPetscObjects();
};
  
}// namespace

#include "data_management/model_order_reduction.tpp"
