#pragma once

#include <vector>

#include <petscmat.h>

class PetscUtility
{
public:
  ///! extract all entries of the PETSc matrix and store values to std::vector
  static void getMatrixEntries(Mat &matrix, std::vector<double> &matrixValues);
  
  ///! extract all entries of the PETSc vector and store values to std::vector
  static void getVectorEntries(Vec &vector, std::vector<double> &vectorValues);
};    