#pragma once

#include <petscmat.h>

#include "control/dihu_context.h"
#include "data_management/data.h"
#include "data_management/model_order_reduction.h"

namespace ModelOrderReduction
{

/** A class for model order reduction techniques.
 */
<<<<<<< HEAD
template<typename FunctionSpaceType>
class MORBase
{
public:
  
  typedef Data::ModelOrderReduction<FunctionSpaceType> Data; //type of Data object
  //! constructor
  MORBase(DihuContext context);
  
  virtual ~MORBase();
  
  //! Set the basis V as Petsc Mat
  void setBasis();
  
  virtual void initialize();
   
protected:
  //! Set the reduced system matrix, A_R=V^T A V
  virtual void setRedSysMatrix(Mat &A, Mat &A_R);
  
  std::shared_ptr<Data> data_;
  
  bool initialized_;
};

}  // namespace

#include "model_order_reduction/mor.tpp"