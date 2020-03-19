#pragma once

#include <petsc.h>
#include <memory>

#include "partition/rank_subset.h"

namespace TimeSteppingScheme
{

namespace NestedMatVecUtility
{

//! from a Petsc Vec with nested type (nestedVec) create a new Petsc Vec (singleVec) that contains all values at once. If the singleVec already exists, do not create again, only copy the values.
void createVecFromNestedVec(Vec nestedVec, Vec &singleVec, std::shared_ptr<Partition::RankSubset> rankSubset);

//! copy the values from a singleVec back to the nested Petsc Vec (nestedVec)
void fillNestedVec(Vec singleVec, Vec nestedVec);

//! from a Petsc Mat with nested type (nestedMat) create a new Petsc Mat (singleMat) that contains all values at once. If the singleMat already exists, do not create again, only copy the values.
void createMatFromNestedMat(Mat nestedMat, Mat &singleMat, std::shared_ptr<Partition::RankSubset> rankSubset);

} 
}  // namespace
