#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include "utility/petsc_utility.h"
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "field_variable/field_variable.h"
#include "partition/partitioned_petsc_vec.h"
#include "partition/partitioned_petsc_mat.h"

namespace Data
{

/**  The datastructures used for timestepping schemes, essentially stores a solution and an increment vector.
 *   nComponents is a constant that determines the number of components of each field variable.
 *   The time stepping works independent of that constant, it is only important for the output files, how many values there are associated with a single node.
 *   E.g. for a cellml model the nComponents should be set to the number of components of this model.
  */
template<typename FunctionSpaceType, int nComponents>
class TimeStepping : public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,nComponents> FieldVariableType;

  //! constructor
  TimeStepping(DihuContext context);

  //! destructur
  ~TimeStepping();

  //! return a reference to the solution vector, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariableType &solution();

  //! return a reference to the increment vector, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariableType &increment();
  
  //! return a reference to the rhs vector, used for the variant 1 of the implicit Euler scheme
  //FieldVariableType &rhs();

  // virtual FieldVariableType &intermediateIncrement() = 0; 
  
  //! get the system matrix required by the implicit time stepping
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> systemMatrix();
  
  //! get the integration matrix if required for the implicit time stepping
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> integrationMatrixRightHandSide();
  
  //! initialize the sytem matrix from a PETSc matrix that was already created, in this case by a MatMatMult
  void initializeSystemMatrix(Mat &systemMatrix);

  //! set the names of the components for the solution field variable
  void setComponentNames(std::vector<std::string> componentNames);
  
  //! perform the final assembly of petsc
  void finalAssembly();
  
  //! print all stored data to stdout
  virtual void print();

  //! return the number of degrees of freedom per mesh node
  static constexpr int getNDofsPerNode();

  //! return the total number of degrees of freedom, this can be a multiple of the number of nodes of the mesh
  virtual dof_no_t nUnknownsLocalWithGhosts();
  
  //! return the total number of degrees of freedom, this can be a multiple of the number of nodes of the mesh
  virtual dof_no_t nUnknownsLocalWithoutGhosts();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,  // geometry
    std::shared_ptr<FieldVariableType>  // solution
  > OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

protected:

  //! initializes the vectors with size
  virtual void createPetscObjects();

  std::shared_ptr<FieldVariableType> solution_;            ///< the vector of the variable of interest
  std::shared_ptr<FieldVariableType> increment_;        ///< the vector for delta u, (note, this might be reduced in future to only a sub-part of the whole data vector if memory consumption is a problem)
  //std::shared_ptr<FieldVariableType> rhs_;     ///for the variant 1 of the implicit Euler scheme
  // std::shared_ptr<FieldVariableType> intermediateIncrement_;
  
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> systemMatrix_;  ///< the system matrix for implicit time stepping
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> integrationMatrixRightHandSide_; ///< the integration matrix for the right hand side if required by the implicit time stepping
  
  std::vector<std::string> componentNames_;      ///< names of the components of the solution and increment variables
  
private:
  //! get maximum number of expected non-zeros in the system matrix
  void getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros);

};

} // namespace Data

#include "data_management/time_stepping.tpp"
