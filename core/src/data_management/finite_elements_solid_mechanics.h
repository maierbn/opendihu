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

class DihuContext;

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
  
  //! return reference to the residual field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,3> &residual();
  
  //! return reference to the increment field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,3> &increment();
  
  //! return reference to the actual geometry field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,3> &geometryActual();
  
  //! return reference to the reference configuration geometry field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,3> &geometryReference();
  
  //! return reference to the displacements field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,3> &displacements();
  
  //! return reference to the rhs field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,3> &externalVirtualEnergy();
  
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
  
  //! get the number of unknows in the solution variable which is 3*nNodes
  virtual dof_no_t nUnknowns();
   
  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>,  // geometryReference
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>,  // actual geometry (stored in mesh)
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>,  // displacements
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>,   // residual
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>   // externalVirtualEnergy
  > OutputFieldVariables;
  
  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();
  
  //! return reference to a stiffness matrix
  Mat &stiffnessMatrix(){LOG(FATAL)<<"this should not be in use";}
  
  //! unneeded getter method
  FieldVariable::FieldVariable<BasisOnMeshType,3> solution(){LOG(FATAL)<<"this should not be in use";}
  FieldVariable::FieldVariable<BasisOnMeshType,3> rightHandSide(){LOG(FATAL)<<"this should not be in use";}
  
private:
 
  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();
 
  //! get maximum number of expected non-zeros in stiffness matrix
  void getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros);

  Mat tangentStiffnessMatrix_;     ///< the tangent stiffness matrix used as jacobian in the nonlinear solver
  Mat massMatrix_;  ///< a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation
  
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> residual_;           ///< the residual vector, needed in the solution process by PETSc
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> increment_;           ///< the increments vector
  
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> geometryReference_;   //< geometry field in reference configuration, the geometry in actual configuration is stored by mesh_
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> displacements_;        //< current displacements
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> externalVirtualEnergy_;        //< the rhs of the nonlinear equation, or δW_ext
  
  bool massMatrixInitialized_ = false;    ///< if the discretization matrix was initialized

};

}  // namespace

#include "data_management/finite_elements_solid_mechanics.tpp"
