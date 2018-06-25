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
Data<BasisOnMeshType>::
Data(DihuContext context) :
  context_(context), rankSubset_(nullptr)
{
}

template<typename BasisOnMeshType>
Data<BasisOnMeshType>::
~Data()
{
}

template<typename BasisOnMeshType>
dof_no_t Data<BasisOnMeshType>::
nUnknowns()
{
  return this->mesh_->nDofs() * 1;  // value for 1 component, can be overloaded to also have the factor nComponents in it
}

template<typename BasisOnMeshType>
void Data<BasisOnMeshType>::
setMesh(std::shared_ptr<BasisOnMeshType> mesh)
{
  this->mesh_ = mesh;
}

template<typename BasisOnMeshType>
void Data<BasisOnMeshType>::
initialize()
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
void Data<BasisOnMeshType>::
setRankSubset(Partition::RankSubset rankSubset)
{
  rankSubset_ = std::make_shared<Partition::RankSubset>(rankSubset);
}
  
template<typename BasisOnMeshType>
const std::shared_ptr<BasisOnMeshType> Data<BasisOnMeshType>::
mesh() const
{
  return this->mesh_;
}

} // namespace Data
