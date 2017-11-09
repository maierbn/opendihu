#pragma once

#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"


class DihuContext;

namespace Data
{

class FiniteElements : public Data
{
public:
 
  //! constructor
  FiniteElements(DihuContext &context);
  
  //! destructur
  ~FiniteElements();
 
  //! return reference to a stiffness matrix
  Mat &stiffnessMatrix();
  
  //! return reference to a right hand side vector
  Vec &rightHandSide();
  
  //! return reference to solution of the system
  Vec &solution();
  
  //! perform the final assembly of petsc
  void finalAssembly();
  
  //! print all stored data to stdout
  void print();
  
  //! if the discretization matrix is already initialized
  bool discretizationMatrixInitialized();
  
  //! create PETSc matrix
  void initializeDiscretizationMatrix();
  
  //! return a reference to the discretization matrix
  Mat &discretizationMatrix();
  
private:
 
  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();
 
  Mat stiffnessMatrix_;     ///< the standard stiffness matrix of the finite element formulation
  Vec rhs_;                 ///< the rhs vector in weak formulation
  Vec solution_;            ///< the vector of the quantity of interest, e.g. displacement
  Mat discretizationMatrix_;  ///< a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation
  
  bool disablePrinting_ = false;    ///< if printing of matrix and vectors is disabled
  bool disableMatrixPrinting_ = false; ///< if the matrix should not be printed
  
  bool discretizationMatrixInitialized_ = false;    ///< if the discretization matrix was initialized
};

}