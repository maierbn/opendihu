#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>

namespace PetscUtility
{

//! extract all entries of the PETSc matrix and store values row-major to std::vector. The vector is resized to the needed size.
void getMatrixEntries(const Mat &matrix, std::vector<double> &matrixValues);

//! extract all entries of the PETSc vector and store values to std::vector
void getVectorEntries(const Vec &vector, std::vector<double> &vectorValues);

//! fills an already existing petsc vector that has the proper size with values
void setVector(const std::vector<double> &vectorValues, Vec &vector);

//! create a Vec Petsc object
void createVector(Vec &vector, int nEntries, std::string name="");

//! return a string of the matrix that can be printed
std::string getStringMatrix(const Mat &matrix);

//! get the string representation of a matrix, from the values
std::string getStringMatrix(std::vector<double> &matrixValues, int nRows, int nColumns, int nRowsGlobal, int nColumnsGlobal);

//! return a string of the matrix and vector, that can be printed
std::string getStringMatrixVector(const Mat &matrix, const Vec &vector);

//! return a string of the vector to be printed
std::string getStringVector(const Vec &vector);

//! return a string representing the sparsity pattern
std::string getStringSparsityPattern(const Mat &matrix);

//! return a string description of the reason why the solution has finished
std::string getStringLinearConvergedReason(KSPConvergedReason convergedReason);

//! return a string description of the reason why the solution has finished
std::string getStringNonlinearConvergedReason(SNESConvergedReason convergedReason);

//! check if the matrix and vector number of entries are correct such that stiffnessMatrix can be multiplied to rhs
void checkDimensionsMatrixVector(Mat &matrix, Vec &input);

}; // namespace

