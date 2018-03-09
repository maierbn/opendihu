#include "data_management/data.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"

namespace Data
{
  
template<typename BasisOnMeshType>
Data<BasisOnMeshType>::Data(DihuContext context) : 
  context_(context), nComponentsPerNode_(1)
{
}

template<typename BasisOnMeshType>
Data<BasisOnMeshType>::~Data()
{
}

template<typename BasisOnMeshType>
void Data<BasisOnMeshType>::setNComponentsPerNode(int n)
{
  this->nComponentsPerNode_ = n;
}

template<typename BasisOnMeshType>
int Data<BasisOnMeshType>::nComponentsPerNode()
{
  return this->nComponentsPerNode_;
}

template<typename BasisOnMeshType>
int Data<BasisOnMeshType>::nDegreesOfFreedom()
{
  return this->mesh_->nNodes() * this->nComponentsPerNode_;
}

template<typename BasisOnMeshType>
void Data<BasisOnMeshType>::setMesh(std::shared_ptr<BasisOnMeshType> mesh)
{
  this->mesh_ = mesh;
}

template<typename BasisOnMeshType>
void Data<BasisOnMeshType>::initialize()
{
  if (!this->initialized_)
  {
    this->createPetscObjects();
    this->initialized_ = true;
  }
  else
  {
    LOG(WARNING) << "Initialize, mesh is already assigned";
  }
}

template<typename BasisOnMeshType>
const std::shared_ptr<BasisOnMeshType> Data<BasisOnMeshType>::mesh() const
{
  return this->mesh_;
}

} // namespace Data
