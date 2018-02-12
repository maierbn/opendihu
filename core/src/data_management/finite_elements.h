#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "basis_on_mesh/mixed_basis_on_mesh.h"

class DihuContext;
//class BasisOnMesh::Mixed;

namespace Data
{

template<typename BasisOnMeshType>
class FiniteElements : public Data<BasisOnMeshType>
{
public:
 
  //! constructor
  FiniteElements(const DihuContext &context);
  
  //! destructor
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
  
  //! get pointers to all field variables that can be written by output writers
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables();
  
private:
 
  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();
 
  //! get maximum number of expected non-zeros in stiffness matrix
  void getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros);

  Mat stiffnessMatrix_;     ///< the standard stiffness matrix of the finite element formulation
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>> rhs_;                 ///< the rhs vector in weak formulation
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>> solution_;            ///< the vector of the quantity of interest, e.g. displacement
  Mat discretizationMatrix_;  ///< a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation
  
  bool disablePrinting_ = false;    ///< if printing of matrix and vectors is disabled
  bool disableMatrixPrinting_ = false; ///< if the matrix should not be printed
  
  bool discretizationMatrixInitialized_ = false;    ///< if the discretization matrix was initialized
};

/** partial specialization for mixed formulation
 */
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
class FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> :
  public Data<HighOrderBasisOnMeshType>
{
public:
 
  //! constructor
  FiniteElements(const DihuContext &context);
  
  //! destructor
  ~FiniteElements();
 
  //! initialize data fields, has to be called after constructor, typically at the begin of a run() method
  void initialize();
  
  //! return reference to a stiffness matrix
  Mat &stiffnessMatrix();
  
  //! return reference to a right hand side vector, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<HighOrderBasisOnMeshType> &rightHandSide();
  
  //! return reference to solution of the system, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<HighOrderBasisOnMeshType> &solution();
  
  //! return reference to geometry in reference configuration field
  FieldVariable::FieldVariable<HighOrderBasisOnMeshType> &geometryReference();
  
  //! return reference to displacement field
  FieldVariable::FieldVariable<HighOrderBasisOnMeshType> &displacement();
  
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
  
  //! get pointers to all field variables that can be written by output writers
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType>>> fieldVariables();
  
  //! set the mixed mesh
  void setMesh(std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> mixedMesh);
  
  //! get the mixed mesh
  std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> mixedMesh();
  
private:
 
  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();
 
  //! get maximum number of expected non-zeros in stiffness matrix
  void getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros);

  //! initialize the geometryReference field variable from the geometry field 
  void initializeFieldVariables();
  
  Mat stiffnessMatrix_;     ///< the standard stiffness matrix of the finite element formulation
  Mat discretizationMatrix_;  ///< a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation
  
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType>> rhs_;                 ///< the rhs vector in weak formulation
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType>> solution_;            ///< the vector of the quantity of interest, e.g. displacement
  
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType>> geometryReference_;   ///< geometry in reference configuration
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType>> displacement_;   ///< displacement fields u
  std::shared_ptr<FieldVariable::FieldVariable<LowOrderBasisOnMeshType>> pressure_;   ///< pressure field p  
  
  
  bool disablePrinting_ = false;    ///< if printing of matrix and vectors is disabled
  bool disableMatrixPrinting_ = false; ///< if the matrix should not be printed
  
  bool discretizationMatrixInitialized_ = false;    ///< if the discretization matrix was initialized
  
  std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> mixedMesh_;   ///< the BasisOnMesh object of Type BasisOnMesh::Mixed
};

}  // namespace

#include "data_management/finite_elements.tpp"
#include "data_management/finite_elements_mixed.tpp"
