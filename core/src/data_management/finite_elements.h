#pragma once

#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"

class DihuContext;

namespace Data
{

template<typename BasisOnMeshType>
class FiniteElements : public Data<BasisOnMeshType>
{
public:
 
  //! constructor
  FiniteElements(const DihuContext &context);
  
  //! destructur
  ~FiniteElements();
 
  //! return reference to a stiffness matrix
  Mat &stiffnessMatrix();
  
  //! return reference to a right hand side vector, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType> &rightHandSide();
  
  //! return reference to solution of the system, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType> &solution();
  
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
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>> rhs_;                 ///< the rhs vector in weak formulation
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>> solution_;            ///< the vector of the quantity of interest, e.g. displacement
  Mat discretizationMatrix_;  ///< a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation
  
  bool disablePrinting_ = false;    ///< if printing of matrix and vectors is disabled
  bool disableMatrixPrinting_ = false; ///< if the matrix should not be printed
  
  bool discretizationMatrixInitialized_ = false;    ///< if the discretization matrix was initialized
};

}  // namespace

#include "data_management/finite_elements.tpp"