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
  
  ///! fills an already existing petsc vector that has the proper size with values
  static void setVector(std::vector<double> &vectorValues, Vec &vector);
  
  ///! create a Vec Petsc object
  static void createVector(Vec &vector, int nEntries, std::string name="");
  
  ///! return a string of the matrix that can be printed
  static std::string getStringMatrix(Mat &matrix);
  
  ///! return a string of the matrix and vector, that can be printed
  static std::string getStringMatrixVector(Mat &matrix, Vec &vector);
  
  ///! return a string of the vector to be printed
  static std::string getStringVector(Vec &vector);
  
  ///! return a string representing the sparsity pattern
  static std::string getStringSparsityPattern(Mat &matrix);
};    