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

template<typename FunctionSpaceType,typename BaseDataType>
StreamlineTracer<FunctionSpaceType,BaseDataType>::
StreamlineTracer(DihuContext context) : Data<FunctionSpaceType>(context), fiberNo_(0)
{
}

template<typename FunctionSpaceType,typename BaseDataType>
StreamlineTracer<FunctionSpaceType,BaseDataType>::
~StreamlineTracer()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType,typename BaseDataType>
void StreamlineTracer<FunctionSpaceType,BaseDataType>::
setBaseData(std::shared_ptr<BaseDataType> baseData)
{
  baseData_ = baseData;

  // set function space
  this->setFunctionSpace(baseData_->functionSpace());
}

template<typename FunctionSpaceType,typename BaseDataType>
void StreamlineTracer<FunctionSpaceType,BaseDataType>::
createPetscObjects()
{
  LOG(DEBUG) << "StreamlineTracer<FunctionSpaceType,BaseDataType>::createPetscObjects()" << std::endl;
  assert(this->functionSpace_);
  
  // create field variables on local partition
  this->gradient_ = this->functionSpace_->template createFieldVariable<3>("gradient");
}

template<typename FunctionSpaceType,typename BaseDataType>
void StreamlineTracer<FunctionSpaceType,BaseDataType>::
createfiberMesh(const std::vector<Vec3> &nodePositions)
{
  std::shared_ptr<FunctionSpacefiber> meshPtr;
 
  // create name for fiber mesh 
  std::stringstream name;
  name << "fiber" << std::setw(5) << std::setfill('0') << fiberNo_;
  fiberNo_++;
  
  // set number of elements. We have 1D linear Lagrange elements, i.e. 2 nodes per element
  const int nElements = nodePositions.size()-1;
  std::array<element_no_t,1> nElementsPerCoordinateDirection{nElements}; 
  
  // create mesh by meshManager
  meshPtr = this->context_.meshManager()->template createFunctionSpace<FunctionSpacefiber>(name.str(), nodePositions, nElementsPerCoordinateDirection);
  
  // get geometry field 
  std::shared_ptr<FieldVariablefiberGeometry> geometryField = std::make_shared<FieldVariablefiberGeometry>(meshPtr->geometryField());
  
  // add geometry field
  this->fiberGeometry_.push_back(geometryField);
}

template<typename FunctionSpaceType,typename BaseDataType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> StreamlineTracer<FunctionSpaceType,BaseDataType>::
gradient()
{
  return this->gradient_;
}

template<typename FunctionSpaceType,typename BaseDataType>
dof_no_t StreamlineTracer<FunctionSpaceType,BaseDataType>::
nNodesLocalWithGhosts()
{
  return this->functionSpace_->nNodesLocalWithGhosts();
}

template<typename FunctionSpaceType,typename BaseDataType>
dof_no_t StreamlineTracer<FunctionSpaceType,BaseDataType>::
nNodesLocalWithoutGhosts()
{
  return this->functionSpace_->nNodesLocalWithoutGhosts();
}

template<typename FunctionSpaceType,typename BaseDataType>
constexpr int StreamlineTracer<FunctionSpaceType,BaseDataType>::
getNDofsPerNode()
{
  return 1;
}

template<typename FunctionSpaceType,typename BaseDataType>
void StreamlineTracer<FunctionSpaceType,BaseDataType>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->gradient_;
  VLOG(4) << "======================";
}

template<typename FunctionSpaceType,typename BaseDataType>
typename StreamlineTracer<FunctionSpaceType,BaseDataType>::OutputFieldVariables StreamlineTracer<FunctionSpaceType,BaseDataType>::
getOutputFieldVariables()
{
  return std::tuple_cat(baseData_->getOutputFieldVariables(),
                        std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>>(gradient_),
                        std::tuple<std::vector<std::shared_ptr<FieldVariablefiberGeometry>>>(fiberGeometry_)
                       );
}


} // namespace
