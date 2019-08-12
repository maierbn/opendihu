#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>

// color codes: https://github.com/shiena/ansicolor/blob/master/README.md
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_LIGHT_GRAY    "\x1b[90m"
#define ANSI_COLOR_LIGHT_WHITE    "\x1b[97m"
#define ANSI_COLOR_RESET   "\x1b[0m"

namespace PetscUtility
{

//! extract all entries of the PETSc matrix and store values row-major to std::vector. The vector is resized to the needed size.
void getMatrixEntries(const Mat &matrix, std::vector<double> &matrixValues);

//! extract all entries of the PETSc vector and store values to std::vector
void getVectorEntries(const Vec &vector, std::vector<double> &vectorValues);

//! fills an already existing petsc vector that has the proper size with values
void setVector(const std::vector<double> &vectorValues, Vec &vector);

//! create a Vec Petsc object, this function is to be removed, currently only needed because of solid mechanics
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

//! write the vector to a file using PetscViewer, format is "default", "ascii" or "matlab"
void dumpVector(std::string filename, std::string format, Vec &vector, MPI_Comm mpiCommunicator, int componentNo=0, int nComponents=1);

//! dump the matrix to a file using PetscViewer, format is "default", "ascii" or "matlab"
//! if the matrix is nested, only one file is created that contains the whole matrix data.
void dumpMatrix(std::string filename, std::string format, Mat &matrix, MPI_Comm mpiCommunicator);

}  // namespace

