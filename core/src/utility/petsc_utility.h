#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>

class PetscUtility
{
public:
  ///! extract all entries of the PETSc matrix and store values row-major to std::vector. The vector is resized to the needed size.
  static void getMatrixEntries(const Mat &matrix, std::vector<double> &matrixValues);
  
  ///! extract all entries of the PETSc vector and store values to std::vector
  static void getVectorEntries(const Vec &vector, std::vector<double> &vectorValues);
  
  ///! fills an already existing petsc vector that has the proper size with values
  static void setVector(const std::vector<double> &vectorValues, Vec &vector);
  
  ///! create a Vec Petsc object
  static void createVector(Vec &vector, int nEntries, std::string name="");
  
  ///! return a string of the matrix that can be printed
  static std::string getStringMatrix(const Mat &matrix);
  
  ///! return a string of the matrix and vector, that can be printed
  static std::string getStringMatrixVector(const Mat &matrix, const Vec &vector);
  
  ///! return a string of the vector to be printed
  static std::string getStringVector(const Vec &vector);
  
  ///! return a string representing the sparsity pattern
  static std::string getStringSparsityPattern(const Mat &matrix);
  
  ///! return a string description of the reason why the solution has finished
  static std::string getStringLinearConvergedReason(KSPConvergedReason convergedReason);
  
  ///! return a string description of the reason why the solution has finished
  static std::string getStringNonlinearConvergedReason(SNESConvergedReason convergedReason);
};    
