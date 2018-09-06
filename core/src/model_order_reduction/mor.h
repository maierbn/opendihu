
#include "control/dihu_context.h"
#include <petscmat.h>

namespace ModelOrderReduction
{

/** A class for model order reduction techniques.
 */
class MORBase
{
public:
  //! constructor
  MORBase(DihuContext context, int n, int k);
  
  //! Basis for the reduced solution, V
  Mat &basis();
  //! Transpose of the basis, V ^T
  Mat &basisTransp();
  //! Reduced system matrix, A_R
  Mat &redSysMatrix();
  
  //! Set the basis V as Petsc Mat
  void setBasis();
   
protected:
  //! Set the reduced system matrix, A_R=V^T A V
  virtual void setRedSysMatrix(Mat A, Mat A_R);
  
private:
  //! Create the matrices and vectors for model order reduction
  void createPetscObjects();
   
  Mat basis_; // V
  Mat basisTransp_; // V^T
  Mat redSysMatrix_;
  
  Vec solution_; // full order solution
  Vec redSolution_; //reduced solution
 
};

}  // namespace