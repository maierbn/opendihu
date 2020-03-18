#include "specialized_solver/multidomain_solver/nested_mat_vec_utility.h"

namespace TimeSteppingScheme
{

namespace NestedMatVecUtility
{

//! from a Petsc Vec with nested type (nestedVec) create a new Petsc Vec (singleVec) that contains all values at once. If the singleVec already exists, do not create again, only copy the values.
void createVecFromNestedVec(Vec nestedVec, Vec &singleVec)
{
  // if Vec object does not yet exist, create new one
  if (singleVec == PETSC_NULL)
  {

  }

  singleVec = nestedVec;
}

//! copy the values from a singleVec back to the nested Petsc Vec (nestedVec)
void fillNestedVec(Vec singleVec, Vec nestedVec)
{

}

//! from a Petsc Mat with nested type (nestedMat) create a new Petsc Mat (singleMat) that contains all values at once. If the singleMat already exists, do not create again, only copy the values.
void createMatFromNestedMat(Mat nestedMat, Mat &singleMat)
{
  // if Mat object does not yet exist, create new one
  if (singleMat == PETSC_NULL)
  {

  }

  singleMat = nestedMat;
}

} 
}  // namespace