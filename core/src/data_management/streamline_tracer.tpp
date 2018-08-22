#include "data_management/time_stepping.h"

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

template<typename BasisOnMeshType,typename BaseDataType>
StreamlineTracer<BasisOnMeshType,BaseDataType>::
StreamlineTracer(DihuContext context) : Data<BasisOnMeshType>(context), fibreNo_(0)
{
}

template<typename BasisOnMeshType,typename BaseDataType>
StreamlineTracer<BasisOnMeshType,BaseDataType>::
~StreamlineTracer()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,typename BaseDataType>
void StreamlineTracer<BasisOnMeshType,BaseDataType>::
setBaseData(std::shared_ptr<BaseDataType> baseData)
{
  baseData_ = baseData;

  // set mesh
  this->setMesh(baseData_->mesh());
}

template<typename BasisOnMeshType,typename BaseDataType>
void StreamlineTracer<BasisOnMeshType,BaseDataType>::
createPetscObjects()
{
  LOG(DEBUG)<<"StreamlineTracer<BasisOnMeshType,BaseDataType>::createPetscObjects()"<<std::endl;
  assert(this->mesh_);
  
  // create partitioning
  Partition::MeshPartition<BasisOnMeshType> partition = this->context_.template createPartitioning<BasisOnMeshType>(this->rankSubset_, this->mesh_);
  
  // create field variables on local partition
  this->gradient_ = this->mesh_->template createFieldVariable<3>("gradient", partition);
}

template<typename BasisOnMeshType,typename BaseDataType>
void StreamlineTracer<BasisOnMeshType,BaseDataType>::
createFibreMesh(const std::vector<Vec3> &nodePositions)
{
  std::shared_ptr<MeshFibre> meshPtr; 
 
  // create name for fibre mesh 
  std::stringstream name;
  name << "Fibre" << std::setw(5) << std::setfill('0') << fibreNo_;
  fibreNo_++;
  
  // set number of elements. We have 1D linear Lagrange elements, i.e. 2 nodes per element
  const int nElements = nodePositions.size()-1;
  std::array<element_no_t,1> nElementsPerCoordinateDirection{nElements}; 
  
  // create mesh by meshManager
  meshPtr = std::static_pointer_cast<MeshFibre>(
     this->context_.meshManager()->template createMesh<MeshFibre>(name.str(), nodePositions, nElementsPerCoordinateDirection));
  
  // get geometry field 
  std::shared_ptr<FieldVariableFibreGeometry> geometryField = std::make_shared<FieldVariableFibreGeometry>(meshPtr->geometryField());
  
  // add geometry field
  this->fibreGeometry_.push_back(geometryField);
}

template<typename BasisOnMeshType,typename BaseDataType>
FieldVariable::FieldVariable<BasisOnMeshType,3> &StreamlineTracer<BasisOnMeshType,BaseDataType>::
gradient()
{
  return *this->gradient_;
}

template<typename BasisOnMeshType,typename BaseDataType>
dof_no_t StreamlineTracer<BasisOnMeshType,BaseDataType>::
nNodesLocalWithGhosts()
{
  return this->mesh_->nNodesLocalWithGhosts();
}

template<typename BasisOnMeshType,typename BaseDataType>
dof_no_t StreamlineTracer<BasisOnMeshType,BaseDataType>::
nNodesLocalWithoutGhosts()
{
  return this->mesh_->nNodesLocalWithoutGhosts();
}

template<typename BasisOnMeshType,typename BaseDataType>
constexpr int StreamlineTracer<BasisOnMeshType,BaseDataType>::
getNDofsPerNode()
{
  return 1;
}

template<typename BasisOnMeshType,typename BaseDataType>
void StreamlineTracer<BasisOnMeshType,BaseDataType>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->gradient_;
  VLOG(4) << "======================";
}

template<typename BasisOnMeshType,typename BaseDataType>
typename StreamlineTracer<BasisOnMeshType,BaseDataType>::OutputFieldVariables StreamlineTracer<BasisOnMeshType,BaseDataType>::
getOutputFieldVariables()
{
  return std::tuple_cat(baseData_->getOutputFieldVariables(),
                        std::tuple<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>>(gradient_),
                        std::tuple<std::vector<std::shared_ptr<FieldVariableFibreGeometry>>>(fibreGeometry_)
                       );
}


} // namespace