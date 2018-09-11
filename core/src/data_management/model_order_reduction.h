#pragma once


namespace Data{

template<typename FunctionSpaceType>  
class ModelOrderReduction:
  public Data<FunctionSpaceType>
{
public:
  //! constructor
  ModelOrderReduction(DihuContext context);
   
  virtual ~ModelOrderReduction();
   
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
   
  //! solution which is of course the reduced solution
  Vec &Solution();
   
  //! The reduced order increment
  Vec &Increment();
   
  virtual void initialize() override;
   
private:
   
  Mat basis_; // V
  Mat basisTransp_; // V^T
  Mat redSysMatrix_;
   
  Vec redSolution_; //reduced solution
  Vec redIncrement_; //reduced increment
   
  //! Create the matrices and vectors for model order reduction
  void createPetscObjects();
};
  
}// namespace

#include "data_management/model_order_reduction.tpp"