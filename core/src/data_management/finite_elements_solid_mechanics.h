#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <tuple>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "equation/type_traits.h"

namespace Data
{
 
template<typename BasisOnMeshType,typename Term>
class FiniteElements<
  BasisOnMeshType,
  Term,
  Equation::isSolidMechanics<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
> : 
  public Data<BasisOnMeshType>
{
public:
 
  //! constructor
  FiniteElements(DihuContext context);
  
  //! destructor
  ~FiniteElements();
 
  //! initialize the object, create all stored data
  virtual void initialize() override;
  
  //! return reference to a stiffness matrix
  Mat &tangentStiffnessMatrix();
  
  //! return reference to the reduced tangent stiffness matrix
  Mat &tangentStiffnessMatrixReduced();
  
  //! return reference to the residual field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &residual();
  
  //! return reference to the increment field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &increment(){LOG(FATAL)<<"this should not be in use";}
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &solution(){LOG(FATAL)<<"this should not be in use";}
  
  //! return reference to the actual geometry field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,3> &geometryActual();
  
  //! return reference to the reference configuration geometry field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,3> &geometryReference();
  
  //! return reference to the displacements field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &displacements();
  
  //! return reference to the wExt field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &externalVirtualWork();
  
  //! return reference to the wInt field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &internalVirtualWork();
  
  //! alias for externalVirtualWork, needed such that rhs setting functionality works
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &rightHandSide();
  
  //! vector to be used for reduced computation, for displacements
  Vec &displacementsReduced();
  
  //! vector to be used for reduced computation, for rhs
  Vec &rhsReduced();
  
  //! vector to be used for reduced computation, for residual in SNES
  Vec &residualReduced();
  
  //! 3D vector to be used in 2D problems for adding to actual geometry which is also 3D
  Vec &fullIncrement();
  
  //! perform the final assembly of petsc
  void finalAssembly();
  
  //! print all stored data to stdout
  void print();
  
  //! if the discretization matrix is already initialized
  bool massMatrixInitialized();
  
  //! create PETSc matrix
  void initializeMassMatrix();
  
  //! if the external virtual energy does not depend on deformation and thus is constant
  bool &externalVirtualWorkIsConstant();
  
  //! return a reference to the discretization matrix
  Mat &massMatrix();
  
  //! get the number of unknows in the solution variable which is 3*nNodes
  virtual dof_no_t nUnknowns();
   
  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>,  // geometryReference
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>,  // actual geometry (stored in mesh)
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>>,  // displacements
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>>,   // residual
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>>   // externalVirtualWork
  > OutputFieldVariables;
  
  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();
  
  //! return reference to a stiffness matrix
  Mat &stiffnessMatrix(){LOG(FATAL)<<"this should not be in use";}

  //! return the value of computeWithReducedVectors. If the vector of unknowns only contains the real degrees of freedom and not the variables with Dirichlet BCs. This is maybe slower because copying of data is required, but the system to solve is smaller
  bool computeWithReducedVectors();
  
  //! initialize the reduced size variables
  void initializeReducedVariables(dof_no_t nDofsReduced);
   
private:
 
  //! initializes the vector and stiffness matrix with size
  void createPetscObjects();
 
  //! get maximum number of expected non-zeros in stiffness matrix
  void getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros);

  Mat tangentStiffnessMatrix_;     ///< the tangent stiffness matrix used as jacobian in the nonlinear solver
  Mat massMatrix_;  ///< a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation
  Mat tangentStiffnessMatrixReduced_;   ///< the tangent stiffness matrix without entries for Dirichlet BC
  
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>> residual_;           ///< the residual vector, needed in the solution process by PETSc
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>> increment_;           ///< the increments vector
  
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> geometryReference_;   //< geometry field in reference configuration, the geometry in actual configuration is stored by mesh_
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>> displacements_;        //< current displacements
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>> externalVirtualWork_;        //< the external virtual work vector δW_ext
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>> internalVirtualWork_;        //< the internal virtual work vector δW_int
  
  Vec displacementsReduced_;   ///< if computeWithReducedVectors_ this vector is used to store the reduced displacements vector that does contain the same is displacements_ but without values for Dirichlet BCs
  Vec rhsReduced_;  ///<  if computeWithReducedVectors_ this vector is used to store the values for the rhs (δW_int - δW_ext) without values that have Dirichlet BCs
  Vec residualReduced_; ///< if computeWithReducedVectors_ this vector is used to store a reduced version of the residual for the nonlinear solver
  Vec fullIncrement_;   ///< only for 2D problems this vec is a 3D vector that is filled from the 2D displacements vector and afterwards added to the geometry values.
  
  bool externalVirtualWorkIsConstant_;   //< if the external virtual energy does not depend on deformation and thus is constant
  
  bool massMatrixInitialized_ = false;    ///< if the discretization matrix was initialized
  const bool computeWithReducedVectors_;    ///< if the vector of unknowns only contains the real degrees of freedom and not the variables with Dirichlet BCs. This is maybe slower because copying of data is required, but the system to solve is smaller.

};

}  // namespace

#include "data_management/finite_elements_solid_mechanics.tpp"
