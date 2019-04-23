#include "field_variable/structured/02b_field_variable_data_structured_for_surface.h"

#include <sstream>
#include <petsc.h>

#include "partition/rank_subset.h"
#include "mesh/mesh_manager.h"

namespace FieldVariable
{

template<typename BasisFunctionType, int nComponents>
FieldVariableDataStructuredForSurface<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::
FieldVariableDataStructuredForSurface(FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents> &rhs, Mesh::face_t face)
{
  // initialize everything from other field variable
  assert(rhs.componentNames_.size() == nComponents);
  std::copy(rhs.componentNames_.begin(), rhs.componentNames_.end(), this->componentNames_.begin());

  this->name_ = rhs.name_;
  this->isGeometryField_ = rhs.isGeometryField_;

  LOG(DEBUG) << "create 2D surface field variable from 3D field variable \"" << rhs.name_ << "\"";

  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType> FunctionSpace3D;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType> FunctionSpace2D;

  std::shared_ptr<Partition::MeshPartition<FunctionSpace3D>> meshPartition3D = rhs.functionSpace()->meshPartition();

  // construct rankSubset

  // collect rank nos
  int rankNoEnd = meshPartition3D->nRanks(0) * meshPartition3D->nRanks(1);
  int rankNoBegin = 0;
  std::vector<int> rankNos(rankNoEnd);

  if (face == Mesh::face_t::face2Plus)
  {
    rankNoBegin = meshPartition3D->nRanks(0) * meshPartition3D->nRanks(1) * (meshPartition3D->nRanks(2) - 1);
  }
  std::iota(rankNos.begin(), rankNos.end(), rankNoBegin);

  std::shared_ptr<Partition::RankSubset> rankSubset = std::make_shared<Partition::RankSubset>(rankNos.begin(), rankNos.end(), rhs.functionSpace->rankSubset());

  LOG(DEBUG) << "rankSubset: " << rankSubset;

  // construct meshPartition
  std::array<node_no_t,2> nElementsLocal{meshPartition3D->nElementsLocal(0), meshPartition3D->nElementsLocal(1)};
  std::array<global_no_t,2> nElementsGlobal{meshPartition3D->nElementsGlobal(0), meshPartition3D->nElementsGlobal(1)};
  std::array<int,2> beginElementGlobal{meshPartition3D->beginElementGlobal_[0], meshPartition3D->beginElementGlobal_[1]};
  std::array<int,2> nRanks{meshPartition3D->nRanks[0], meshPartition3D->nRanks[1]};

  std::shared_ptr<Partition::MeshPartition<FunctionSpace2D>> meshPartition = std::make_shared<Partition::MeshPartition<FunctionSpace2D>>(
    nElementsLocal, nElementsGlobal, beginElementGlobal, nRanks, rankSubset
  );

  LOG(DEBUG) << "meshPartition 3D: " << *meshPartition3D;
  LOG(DEBUG) << "meshPartition 2D: " << *meshPartition;

  // create surface function space
  std::stringstream functionSpaceName;
  functionSpaceName << rhs->functionSpace()->meshName() << "_surface";

  const int nDofsLocalWithoutGhosts = meshPartition->nDofsLocalWithoutGhosts();

  // create node positions (without ghost nodes)
  std::vector<Vec3> localNodePositions2D;

  std::vector<dof_no_t> surfaceDofs(nDofsLocalWithoutGhosts);

  dof_no_t surfaceDofsBegin = 0;
  if (face == Mesh::face_t::face2Plus)
  {
    surfaceDofsBegin = meshPartition3D->nDofsLocalWithoutGhosts() - nDofsLocalWithoutGhosts;
  }

  std::iota(surfaceDofs.begin(), surfaceDofs.end(), surfaceDofsBegin);

  LOG(DEBUG) << "nDofsLocalWithoutGhosts: " << nDofsLocalWithoutGhosts << ", surfaceDofs: " << surfaceDofs;

  rhs.functionSpace()->geometryField()->getValues(surfaceDofs, localNodePositions2D);
  std::array<element_no_t,2> nElementsPerCoordinateDirectionLocal();

  LOG(DEBUG) << "localNodePositions2D: " << localNodePositions2D << ", now create function space";

  this->functionSpace_ = DihuContext::meshManager()->createFunctionSpaceWithGivenMeshPartition<FunctionSpace2D>(
    functionSpaceName.str(), meshPartition, localNodePositions2D, nElementsLocal, nRanks);

  assert(this->functionSpace_);
  assert(this->functionSpace_->meshPartition());

  VLOG(1) << "construct field variable \"" << this->name_ << "\" from other field variable \"" << rhs.name() << "\".";

  // create new distributed petsc vec as copy of rhs values vector
  if (rhs.partitionedPetscVec())
  {
    // if rhs is not a geometry field an therefore has a partitionedPetscVec, use that
    this->values_ = std::make_shared<PartitionedPetscVec<FunctionSpace2D,nComponents>>(*rhs.partitionedPetscVec(), this->name_);
  }
  else
  {
    // else create new values vector
    this->values_ = std::make_shared<PartitionedPetscVec<FunctionSpace2D,nComponents>>(this->functionSpace_->meshPartition(), this->name_);
  }

  LOG(DEBUG) << "copy values from rhs field variable";

  // copy values from rhs
  std::vector<double> values;
  std::vector<int> dofs(nDofsLocalWithoutGhosts);
  std::iota(dofs.begin(), dofs.end(), 0);

  // for every component, get values from rhs and set in values
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    values.clear();
    rhs->getValues(componentNo, nDofsLocalWithoutGhosts, surfaceDofs.data(), values);
    LOG(DEBUG) << "component " << componentNo << ", values: " << values;
    this->values_->setValues(componentNo, nDofsLocalWithoutGhosts, dofs.data(), values, INSERT_VALUES);
  }
}

}  // namespace
