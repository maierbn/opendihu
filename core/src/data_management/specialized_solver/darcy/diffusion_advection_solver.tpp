#include "data_management/specialized_solver/darcy/diffusion_advection_solver.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data
{

template<typename FunctionSpaceType>
DiffusionAdvectionSolver<FunctionSpaceType>::
DiffusionAdvectionSolver(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType>
void DiffusionAdvectionSolver<FunctionSpaceType>::
initialize()
{
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();

  // create th slot connector data object
  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();

  // add all needed field variables to be transferred

  // add component 0 of fieldvariableA_
  slotConnectorData_->addFieldVariable(this->solution_, 0);

  // There is addFieldVariable(...) and addFieldVariable2(...) for the two different field variable types,
  // Refer to "slot_connection/slot_connector_data.h" for details.

  // you can also access settings from the python config here:
  std::string option1 = this->context_.getPythonConfig().getOptionString("option1", "default string");
  LOG(DEBUG) << "In data object, parsed option1: [" << option1 << "].";

  // parse slot names for all slot connector data slots. here we only expect one value because we have set one slot (fieldVariableA)
  this->context_.getPythonConfig().getOptionVector("slotNames", slotConnectorData_->slotNames);

  // make sure that there are as many slot names as slots
  slotConnectorData_->slotNames.resize(slotConnectorData_->nSlots());
}

template<typename FunctionSpaceType>
void DiffusionAdvectionSolver<FunctionSpaceType>::
getPetscMemoryParameters(int &nNonZerosDiagonal, int &nNonZerosOffdiagonal)
{
  const int D = FunctionSpaceType::dim();
  const int nDofsPerNode = FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpaceType::BasisFunction>::nDofsPerNode();
  const int nDofsPerElement = FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpaceType::BasisFunction>::nDofsPerElement();
  const int nOverlaps = (nDofsPerElement*2 - 1) * nDofsPerNode;   // number of nodes of 2 neighbouring 1D elements (=number of ansatz functions in support of center ansatz function)

  // due to PETSc storage nNonZerosDiagonal and nNonZerosOffdiagonal should be both set to the maximum number of non-zero entries per row

  switch (D)
  {
  case 1:
    nNonZerosDiagonal = nOverlaps;
    nNonZerosOffdiagonal = nNonZerosDiagonal;
    break;
  case 2:
    nNonZerosDiagonal = pow(nOverlaps, 2);   // because of boundary conditions there can be more entries, which are all zero, but stored as non-zero
    nNonZerosOffdiagonal = nNonZerosDiagonal;
    break;
  case 3:
    nNonZerosDiagonal = pow(nOverlaps, 3);
    nNonZerosOffdiagonal = nNonZerosDiagonal;
    break;
  };

  LOG(DEBUG) << "nOverlaps = (" << nDofsPerElement << "*2 - 1) * " << nDofsPerNode << " = " << nOverlaps
    << ", nNonZerosDiagonal = " << nOverlaps << "^" << D << " = " << nNonZerosDiagonal << ", 2*" << nNonZerosDiagonal << "=" << 2*nNonZerosDiagonal;
}

template<typename FunctionSpaceType>
void DiffusionAdvectionSolver<FunctionSpaceType>::
createPetscObjects()
{
  assert(this->functionSpace_);

  // get the partitioning from the function space
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition = this->functionSpace_->meshPartition();

  // Here, the actual field variables will be created.
  // Make sure, the number of components matches. The string is the name of the field variable. It will also be used in the VTK output files.
  this->solution_ = this->functionSpace_->template createFieldVariable<1>("solution");
  this->increment_ = this->functionSpace_->template createFieldVariable<1>("increment");

  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int nNonZerosDiagonal = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int nNonZerosOffdiagonal = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(nNonZerosDiagonal, nNonZerosOffdiagonal);

  LOG(DEBUG) << "d=" << this->functionSpace_->dimension()
    << ", number of diagonal non-zeros: " << nNonZerosDiagonal << ", number of off-diagonal non-zeros: " <<nNonZerosOffdiagonal;

  LOG(DEBUG) << "create new stiffnessMatrix";
  this->vMatrix_ = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(meshPartition, 1, nNonZerosDiagonal, nNonZerosOffdiagonal, "vMatrix");
}

// ... add a "getter" method for each fieldvariable with the same name as the field variable (but without underscore)
template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> DiffusionAdvectionSolver<FunctionSpaceType>::
solution()
{
  return this->solution_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> DiffusionAdvectionSolver<FunctionSpaceType>::
increment()
{
  return this->increment_;
}

template<typename FunctionSpaceType>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> DiffusionAdvectionSolver<FunctionSpaceType>::
vMatrix()
{
  return this->vMatrix_;
}

template<typename FunctionSpaceType>
std::shared_ptr<typename DiffusionAdvectionSolver<FunctionSpaceType>::SlotConnectorDataType> DiffusionAdvectionSolver<FunctionSpaceType>::
getSlotConnectorData()
{
  // return the slot connector data object
  return this->slotConnectorData_;
}

template<typename FunctionSpaceType>
typename DiffusionAdvectionSolver<FunctionSpaceType>::FieldVariablesForOutputWriter DiffusionAdvectionSolver<FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files by the output writer

  // get the geometry field, which is always needed, from the function space
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  return std::make_tuple(
    geometryField,
    this->solution_,
    this->increment_   // add all field variables that should appear in the output file. Of course, this list has to match the type in the header file.
  );
}

} // namespace
