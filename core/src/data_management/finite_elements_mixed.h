#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "data_management/finite_elements_mixed.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "basis_on_mesh/mixed_basis_on_mesh.h"

class DihuContext;
//class BasisOnMesh::Mixed;

namespace Data
{
/** partial specialization for mixed formulation
 */
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
class FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term> :
  public Data<HighOrderBasisOnMeshType>,
  public DiffusionTensor<HighOrderBasisOnMeshType::dim()>
{
public:
 
  typedef FieldVariable::FieldVariable<HighOrderBasisOnMeshType,HighOrderBasisOnMeshType::dim()> HighOrderFieldVariableType;
  typedef FieldVariable::FieldVariable<LowOrderBasisOnMeshType,LowOrderBasisOnMeshType::dim()> LowOrderFieldVariableType;
 
  //! constructor
  FiniteElements(DihuContext context);
  
  //! destructor
  ~FiniteElements();
 
  //! initialize data fields, has to be called after constructor, typically at the begin of a run() method
  virtual void initialize() override;
  
  //! return reference to a right hand side vector, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<HighOrderBasisOnMeshType,HighOrderBasisOnMeshType::dim()> &rightHandSide();
  
  //! return reference to solution of the system, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<HighOrderBasisOnMeshType,HighOrderBasisOnMeshType::dim()> &solution();
  
  //! return reference to geometry in reference configuration field
  FieldVariable::FieldVariable<HighOrderBasisOnMeshType,HighOrderBasisOnMeshType::dim()> &geometryReference();
  
  //! return reference to displacement field
  FieldVariable::FieldVariable<HighOrderBasisOnMeshType,HighOrderBasisOnMeshType::dim()> &displacement();
  
  //! return reference to pressure field
  FieldVariable::FieldVariable<LowOrderBasisOnMeshType,1> &pressure();
  
  //! return reference to fu, right hand side in weak form of displacements
  FieldVariable::FieldVariable<HighOrderBasisOnMeshType,1> &f();
  
  //! perform the final assembly of petsc
  void finalAssembly();
  
  //! print all stored data to stdout
  void print();
  
  //! if the discretization matrix is already initialized
  bool massMatrixInitialized();
  
  //! create PETSc matrix
  void initializeMassMatrix();
  
  //! return reference to a stiffness matrix
  Mat &stiffnessMatrix();
  
  //! return a reference to the discretization matrix
  Mat &massMatrix();
  
  //! return a reference to the element stiffness matrix kuu
  Mat &kuu();
  
  //! return a reference to the element stiffness matrix kup
  Mat &kup();
  
  //! return a reference to the element stiffness matrix kpp
  Mat &kpp();

  //! return a reference to the inverse kpp
  Mat &kppInverse();
  
  //! return a reference to the transpose of kup
  Mat &kupTranspose();
  
  //! return a reference to a temporary matrix
  Mat &tempKupMatrix();
  
  //! return a reference to the schur complement kuu-kup*kpp^-1*kup^T
  Mat &schurComplement();
      
  //! get references to fu, fp, tempKppFp, tempKupKppFp vectors
  void getLocalVectors(Vec &fu, Vec &fp, Vec &tempKppFp, Vec &tempKupKppFp);
    
  //! get pointers to all field variables that can be written by output writers
  
  //! set the mixed mesh
  void setMesh(std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> mixedMesh);
  
  //! get the mixed mesh
  std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> mixedMesh();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<HighOrderFieldVariableType>,  // geometry
    std::shared_ptr<HighOrderFieldVariableType>,  // solution
    std::shared_ptr<HighOrderFieldVariableType>   // rhs
  > OutputFieldVariables;
  
  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();
  
private:
 
  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();
 
  //! get maximum number of expected non-zeros in stiffness matrix
  void getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros);

  //! initialize the geometryReference field variable from the geometry field 
  void initializeFieldVariables();
  
  //! initialize the element stiffness matrices kuu, kup, kpp
  void initializeMatrices();
  
  // global matrices
  Mat stiffnessMatrix_;     ///< the standard stiffness matrix of the finite element formulation
  Mat massMatrix_;  ///< a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation
  
  // element-local matrices
  Mat kuu_;   ///< element stiffness matrix d2/(du*du)int Psi
  Mat kup_;   ///< element stiffness matrix d2/(du*dp)int Psi
  Mat kpp_;   ///< element stiffness matrix d2/(dp*dp)int Psi
  Mat kppInverse_;  ///< inverse of kpp
  Mat kupTranspose_; ///< transpose of kup
  Mat tempKupMatrix_; ///< temporary helper matrix
  Mat schurComplement_;  ///< schur complement kuu-kup*kpp^-1*kup^T
  Vec fu_;    ///< element temporary vector for rhs in weak form
  Vec fp_;    ///< element temporary vector for rhs in weak form
  Vec tempKppFp_;  ///< temporary vector
  Vec tempKupKppFp_;  ///< temporary vector
   
  // field variables
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,3>> rhs_;                 ///< the rhs vector in weak formulation
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,3>> solution_;            ///< the vector of the quantity of interest, e.g. displacement
  
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,3>> geometryReference_;   ///< geometry in reference configuration
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,3>> displacement_;   ///< displacement fields u
  std::shared_ptr<FieldVariable::FieldVariable<LowOrderBasisOnMeshType,1>> pressure_;   ///< pressure field p  
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,1>> f_;   ///< right hand side in weak form of displacements (temporary)
  
  bool massMatrixInitialized_ = false;    ///< if the discretization matrix was initialized
  
  std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> mixedMesh_;   ///< the BasisOnMesh object of Type BasisOnMesh::Mixed
};

}  // namespace

#include "data_management/finite_elements_mixed.tpp"
