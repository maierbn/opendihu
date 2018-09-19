#pragma once

#include "function_space/function_space.h"
#include "partition/partitioned_petsc_mat.h"


namespace Data{

template<typename FullFunctionSpaceType>  
class ModelOrderReduction:
  public Data<FunctionSpace::Generic>
{
public:
  typedef FieldVariable::FieldVariable<FunctionSpace::Generic,1> FieldVariableType;
  
  //! constructor
  ModelOrderReduction(DihuContext context);
   
  virtual ~ModelOrderReduction();
  
  //! initialize the full order function space
  virtual void setFullFunctionSpace(std::shared_ptr<FullFunctionSpaceType> mesh);
   
  //! Basis for the reduced solution, V
  Mat &basis();
   
  //! Transpose of the basis, V ^T
  Mat &basisTransp(); 
   
  //! initializes the basis V from an already existant Petsc Mat !?
  //void initializeBasis(Mat &basis);
   
  //! Reduced system matrix, A_R
  Mat &redSysMatrix();
   
  //! initializes the redSysMatrix from an already existant Petsc Mat !?
  //void initializeRedSysMatrix(Mat &A_R);
   
  //! the reduced solution
  std::shared_ptr<FieldVariableType> &redSolution();
   
  //! The reduced order increment
  std::shared_ptr<FieldVariableType> &redIncrement();
   
  virtual void initialize() override;
   
private:
   
  Mat basis_; // V
  Mat basisTransp_; // V^T
  Mat redSysMatrix_;
   
  Vec redSolution_; //reduced solution
  Vec redIncrement_; //reduced increment
  
  std::shared_ptr<FullFunctionSpaceType> fullFunctionspace_;
  
  //! Create the matrices and vectors for model order reduction
  void createPetscObjects();
};
  
}// namespace

#include "data_management/model_order_reduction.tpp"