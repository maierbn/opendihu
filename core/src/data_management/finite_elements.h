#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "data_management/diffusion_tensor.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "basis_on_mesh/mixed_basis_on_mesh.h"

class DihuContext;
//class BasisOnMesh::Mixed;

namespace Data
{
 
template<typename BasisOnMeshType>
class FiniteElements : 
  public Data<BasisOnMeshType>,
  public DiffusionTensor<BasisOnMeshType::dim()>
{
public:
 
  //! constructor
  FiniteElements(DihuContext context);
  
  //! destructor
  ~FiniteElements();
 
  //! initialize the object, create all stored data
  virtual void initialize() override;
  
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
  bool massMatrixInitialized();
  
  //! create PETSc matrix
  void initializeMassMatrix();
  
  //! return a reference to the discretization matrix
  Mat &massMatrix();
  
  //! get pointers to all field variables that can be written by output writers
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables();
  
  void debug(std::string name);
private:
 
  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();
 
  //! get maximum number of expected non-zeros in stiffness matrix
  void getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros);

  Mat stiffnessMatrix_;     ///< the standard stiffness matrix of the finite element formulation
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>> rhs_;                 ///< the rhs vector in weak formulation
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>> solution_;            ///< the vector of the quantity of interest, e.g. displacement
  Mat massMatrix_;  ///< a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation
  
  bool disablePrinting_ = false;    ///< if printing of matrix and vectors is disabled
  bool disableMatrixPrinting_ = false; ///< if the matrix should not be printed
  
  bool massMatrixInitialized_ = false;    ///< if the discretization matrix was initialized
};

}  // namespace

#include "data_management/finite_elements.tpp"
#include "data_management/finite_elements_mixed.h"
