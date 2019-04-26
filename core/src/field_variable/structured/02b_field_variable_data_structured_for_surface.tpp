#include "field_variable/structured/02b_field_variable_data_structured_for_surface.h"

#include <sstream>
#include <petsc.h>

#include "partition/rank_subset.h"
#include "mesh/mesh_manager.h"

namespace FieldVariable
{

template<typename BasisFunctionType, int nComponents>
FieldVariableDataStructuredForSurface<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::
FieldVariableDataStructuredForSurface(FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents> &rhs, Mesh::face_t face, bool &ownRankInvolvedInOutput)
{
  // initialize everything from other field variable, which is 3D but this field variable is 2D
  assert(rhs.getNComponents() == nComponents);
  std::copy(rhs.componentNames().begin(), rhs.componentNames().end(), this->componentNames_.begin());

  this->name_ = rhs.name();
  this->isGeometryField_ = rhs.isGeometryField();

  LOG(DEBUG) << "create 2D surface field variable from 3D field variable \"" << rhs.name() << "\"";

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

  std::shared_ptr<Partition::RankSubset> rankSubset = std::make_shared<Partition::RankSubset>(rankNos.begin(), rankNos.end(), rhs.functionSpace()->meshPartition()->rankSubset());

  LOG(DEBUG) << "created rank subset from ranks " << rankNos << ", own is contained: " << rankSubset->ownRankIsContained();

  ownRankInvolvedInOutput = true;
  if (!rankSubset->ownRankIsContained())
  {
    this->values_ = nullptr;
    ownRankInvolvedInOutput = false;
    return;
  }

  LOG(DEBUG) << "rankSubset: " << *rankSubset;

  // construct meshPartition
  std::array<node_no_t,2> nElementsLocal{meshPartition3D->nElementsLocal(0), meshPartition3D->nElementsLocal(1)};
  std::array<global_no_t,2> nElementsGlobal{(global_no_t)meshPartition3D->nElementsGlobal(0), (global_no_t)meshPartition3D->nElementsGlobal(1)};
  std::array<int,2> beginElementGlobal{meshPartition3D->beginElementGlobal(0), meshPartition3D->beginElementGlobal(1)};
  std::array<int,2> nRanks{meshPartition3D->nRanks(0), meshPartition3D->nRanks(1)};

  std::shared_ptr<Partition::MeshPartition<FunctionSpace2D>> meshPartition = std::make_shared<Partition::MeshPartition<FunctionSpace2D>>(
    nElementsLocal, nElementsGlobal, beginElementGlobal, nRanks, rankSubset
  );

  LOG(DEBUG) << "meshPartition 3D: " << *meshPartition3D;
  LOG(DEBUG) << "meshPartition 2D: " << *meshPartition;
  LOG(DEBUG) << "nElementsGlobal: " << nElementsGlobal;

  // create surface function space
  std::stringstream functionSpaceName;
  functionSpaceName << rhs.functionSpace()->meshName() << "_surface";

  const int nDofsLocalWithoutGhosts = meshPartition->nDofsLocalWithoutGhosts();

  // create node positions (without ghost nodes)
  std::vector<Vec3> localNodePositions2D;


  surfaceDofs_.resize(nDofsLocalWithoutGhosts);

  dof_no_t surfaceDofsBegin = 0;
  if (face == Mesh::face_t::face2Plus)
  {
    surfaceDofsBegin = meshPartition3D->nDofsLocalWithoutGhosts() - nDofsLocalWithoutGhosts;
  }

  std::iota(surfaceDofs_.begin(), surfaceDofs_.end(), surfaceDofsBegin);

  LOG(DEBUG) << "nDofsLocalWithoutGhosts: " << nDofsLocalWithoutGhosts << ", surfaceDofs: " << surfaceDofs_;

  rhs.functionSpace()->geometryField().getValues(surfaceDofs_, localNodePositions2D);
  std::array<element_no_t,2> nElementsPerCoordinateDirectionLocal();

  LOG(DEBUG) << "localNodePositions2D: " << localNodePositions2D << ", now create function space";

  // if function space of type FunctionSpace2D and mesh name "<functionspace3d_name>_surface" already exists (because it was already created, e.g. for a previous output field variable)
  if (DihuContext::meshManager()->hasFunctionSpaceOfType<FunctionSpace2D>(functionSpaceName.str()))
  {
    this->functionSpace_ = DihuContext::meshManager()->functionSpace<FunctionSpace2D>(functionSpaceName.str());
  }
  else
  {
    this->functionSpace_ = DihuContext::meshManager()->createFunctionSpaceWithGivenMeshPartition<FunctionSpace2D>(
      functionSpaceName.str(), meshPartition, localNodePositions2D, nElementsLocal, nRanks);
  }

  assert(this->functionSpace_);
  assert(this->functionSpace_->meshPartition());

  VLOG(1) << "construct field variable \"" << this->name_ << "\" from other field variable \"" << rhs.name() << "\".";

  // create new distributed petsc vec as copy of rhs values vector
  // create new values vector
  this->values_ = std::make_shared<PartitionedPetscVec<FunctionSpace2D,nComponents>>(this->functionSpace_->meshPartition(), this->name_);

  this->setValues(rhs);
}

template<typename BasisFunctionType, int nComponents>
void FieldVariableDataStructuredForSurface<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::
setValues(FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents> &rhs)
{
  LOG(DEBUG) << "copy values from rhs field variable";

  // copy values from rhs
  static std::vector<double> values;

  // for every component, get values from rhs and set in values
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    values.clear();
    rhs.getValues(componentNo, surfaceDofs_, values);

    LOG(DEBUG) << "component " << componentNo << ", values: " << values;
    this->values_->setValues(componentNo, this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts(), this->functionSpace_->meshPartition()->dofNosLocal().data(), values.data(), INSERT_VALUES);
  }

  //this->values_->startGhostManipulation();
  this->values_->finishGhostManipulation();
}

}  // namespace
