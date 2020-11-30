#include "field_variable/structured/02b_field_variable_data_structured_for_surface.h"

#include <sstream>
#include <petsc.h>

#include "partition/rank_subset.h"
#include "mesh/mesh_manager/mesh_manager.h"
#include "mesh/mapping_between_meshes/mapping/02_composite.h"

namespace FieldVariable
{

template<typename BasisFunctionType, int nComponents>
FieldVariableDataStructuredForSurface<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::
FieldVariableDataStructuredForSurface(FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents> &rhs, Mesh::face_t face, bool &ownRankInvolvedInOutput) :
  FieldVariableDataStructured<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::FieldVariableDataStructured()
{
  // initialize everything from other field variable, which is 3D but this field variable is 2D
  assert(rhs.getNComponents() == nComponents);
  std::copy(rhs.componentNames().begin(), rhs.componentNames().end(), this->componentNames_.begin());

  this->name_ = rhs.name();
  this->isGeometryField_ = rhs.isGeometryField();

  LOG(DEBUG) << "create 2D surface field variable from 3D field variable \"" << rhs.name() << "\", face: " << Mesh::getString(face) << ", is geometry: " << this->isGeometryField_;

  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType> FunctionSpace3D;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType> FunctionSpace2D;

  std::shared_ptr<Partition::MeshPartition<FunctionSpace3D>> meshPartition3D = rhs.functionSpace()->meshPartition();

  // construct rankSubset
  std::vector<PetscInt> rankNos;
  std::array<PetscInt,3> nRanksPerCoordinateDirection({meshPartition3D->nRanks(0), meshPartition3D->nRanks(1), meshPartition3D->nRanks(2)});
  getSurfaceNumbers(nRanksPerCoordinateDirection, 1, face, rankNos);

  std::shared_ptr<Partition::RankSubset> rankSubset = std::make_shared<Partition::RankSubset>(rankNos.begin(), rankNos.end(), rhs.functionSpace()->meshPartition()->rankSubset());

  VLOG(1) << "created rank subset from ranks " << rankNos << " " << nRanksPerCoordinateDirection << ", own is contained: " << rankSubset->ownRankIsContained();

  ownRankInvolvedInOutput = true;
  if (!rankSubset->ownRankIsContained())
  {
    this->values_ = nullptr;
    ownRankInvolvedInOutput = false;
    return;
  }

  VLOG(1) << "rankSubset: " << *rankSubset;

  // construct meshPartition
  int dimensionIndex0 = 0;
  int dimensionIndex1 = 1;

  if (face == Mesh::face_t::face0Minus || face == Mesh::face_t::face0Plus)
  {
    dimensionIndex0 = 1;
    dimensionIndex1 = 2;
  }
  else if (face == Mesh::face_t::face1Minus || face == Mesh::face_t::face1Plus)
  {
    dimensionIndex0 = 0;
    dimensionIndex1 = 2;
  }

  std::array<node_no_t,2> nElementsLocal{meshPartition3D->nElementsLocal(dimensionIndex0), meshPartition3D->nElementsLocal(dimensionIndex1)};
  std::array<global_no_t,2> nElementsGlobal{(global_no_t)meshPartition3D->nElementsGlobal(dimensionIndex0), (global_no_t)meshPartition3D->nElementsGlobal(dimensionIndex1)};
  std::array<global_no_t,2> beginElementGlobal{meshPartition3D->beginElementGlobal(dimensionIndex0), meshPartition3D->beginElementGlobal(dimensionIndex1)};
  std::array<int,2> nRanks{meshPartition3D->nRanks(dimensionIndex0), meshPartition3D->nRanks(dimensionIndex1)};

  std::shared_ptr<Partition::MeshPartition<FunctionSpace2D>> meshPartition = std::make_shared<Partition::MeshPartition<FunctionSpace2D>>(
    nElementsLocal, nElementsGlobal, beginElementGlobal, nRanks, rankSubset
  );

  VLOG(1) << "meshPartition 3D: " << *meshPartition3D;
  VLOG(1) << "meshPartition 2D: " << *meshPartition;
  VLOG(1) << "nElementsGlobal: " << nElementsGlobal;

  // create surface function space
  std::stringstream functionSpaceName;
  functionSpaceName << rhs.functionSpace()->meshName() << "_surface" << Mesh::getString(face);

  const dof_no_t nDofsLocalWithoutGhosts = meshPartition->nDofsLocalWithoutGhosts();

  // create node positions (without ghost nodes)
  std::vector<Vec3> localNodePositions2D;

  std::array<PetscInt,3> nNodesPerCoordinateDirection({meshPartition3D->nNodesLocalWithoutGhosts(0), meshPartition3D->nNodesLocalWithoutGhosts(1), meshPartition3D->nNodesLocalWithoutGhosts(2)});
  getSurfaceNumbers(nNodesPerCoordinateDirection, rhs.functionSpace()->nDofsPerNode(), face, surfaceDofs_);

  VLOG(1) << "nDofsLocalWithoutGhosts: " << nDofsLocalWithoutGhosts << ", surfaceDofs: " << surfaceDofs_;

  assert(nDofsLocalWithoutGhosts == surfaceDofs_.size());

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
FieldVariableDataStructuredForSurface<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::
FieldVariableDataStructuredForSurface(FieldVariable<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>,BasisFunctionType>,nComponents> &rhs, Mesh::face_t face, bool &ownRankInvolvedInOutput) :
  FieldVariableDataStructuredForSurface<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::
  FieldVariableDataStructuredForSurface(*rhs.subFieldVariable(-1), face, ownRankInvolvedInOutput)
{
}

template<typename BasisFunctionType, int nComponents>
void FieldVariableDataStructuredForSurface<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::
setValues(FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents> &rhs)
{
  VLOG(1) << "copy values from rhs field variable";

  // if this->values_ is nullptr, this means the current rank is not part of the computation of this face
  if (!this->values_)
    return;

  assert(this->values_);

  // copy values from rhs
  static std::vector<double> values;

  // for every component, get values from rhs and set in values
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    values.clear();
    rhs.getValues(componentNo, surfaceDofs_, values);

    //VLOG(1) << "component " << componentNo << ", values: " << values;
    this->values_->setValues(componentNo, this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts(), this->functionSpace_->meshPartition()->dofNosLocal().data(), values.data(), INSERT_VALUES);
  }

  this->values_->setRepresentationGlobal();
  this->values_->startGhostManipulation();
}

template<typename BasisFunctionType, int nComponents>
void FieldVariableDataStructuredForSurface<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::
setValues(FieldVariable<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>,BasisFunctionType>,nComponents> &rhs)
{
  VLOG(1) << "copy values from rhs field variable (composite)";
  rhs.updateSubFieldVariables();

  // if this->values_ is nullptr, this means the current rank is not part of the computation of this face
  if (!this->values_)
    return;

  assert(this->values_);

  // copy values from rhs
  static std::vector<double> values;

  // for every component, get values from rhs and set in values
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    values.clear();
    // get values from last sub field variable of composite mesh
    rhs.subFieldVariableWithoutUpdate(-1)->getValues(componentNo, surfaceDofs_, values);

    //VLOG(1) << "component " << componentNo << ", values: " << values;
    this->values_->setValues(componentNo, this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts(), this->functionSpace_->meshPartition()->dofNosLocal().data(), values.data(), INSERT_VALUES);
  }

  this->values_->setRepresentationGlobal();
  this->values_->startGhostManipulation();
}

template<typename BasisFunctionType, int nComponents>
void FieldVariableDataStructuredForSurface<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::
getSurfaceNumbers(const std::array<node_no_t,3> size, int nDofsPerNode, Mesh::face_t face, std::vector<node_no_t> &surfaceNumbers)
{

  node_no_t nNumbersTotal;

  switch(face)
  {
  case Mesh::face_t::face0Minus:   // left
  case Mesh::face_t::face0Plus:    // right
    nNumbersTotal = size[1]*size[2];
    break;
  case Mesh::face_t::face1Minus:   // front
  case Mesh::face_t::face1Plus:    // back
    nNumbersTotal = size[0]*size[2];
    break;
  case Mesh::face_t::face2Minus:   // bottom
  case Mesh::face_t::face2Plus:    // top
    nNumbersTotal = size[0]*size[1];
    break;
  default:
    break;
  }
  surfaceNumbers.resize(nNumbersTotal*nDofsPerNode);

  VLOG(1) << "size: " << size << ", nNumbersTotal: " << nNumbersTotal << ", nDofsPerNode: " << nDofsPerNode;

  // collect rank nos
  if (face == Mesh::face_t::face1Minus)   // front
  {
    int i = 0;
    int rankNo = 0;
    const dof_no_t rankNoStride = size[0] * (size[1] - 1);
    for (int zIndex = 0; zIndex < size[2]; zIndex++, rankNo += rankNoStride)
    {
      for (int xIndex = 0; xIndex < size[0]; xIndex++, rankNo++)
      {
        for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++, i++)
        {
          surfaceNumbers[i] = rankNo * nDofsPerNode + nodalDofIndex;
        }
      }
    }
  }
  else if (face == Mesh::face_t::face1Plus)   // back
  {
    int i = 0;
    const int rankNoStride = size[0] * (size[1] - 1);
    int rankNo = rankNoStride;
    for (int zIndex = 0; zIndex < size[2]; zIndex++, rankNo += rankNoStride)
    {
      for (int xIndex = 0; xIndex < size[0]; xIndex++, rankNo++)
      {
        for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++, i++)
        {
          surfaceNumbers[i] = rankNo * nDofsPerNode + nodalDofIndex;
        }
      }
    }
  }
  else
  {
    node_no_t rankNoEnd = size[0]*size[1]*size[2];
    int rankNoBegin = 0;
    int rankNoStride = 1;

    switch(face)
    {
    case Mesh::face_t::face0Minus:   // left
      rankNoStride = size[0];
      break;
    case Mesh::face_t::face0Plus:   // right
      rankNoBegin = size[0]-1;
      rankNoStride = size[0];
      break;
    case Mesh::face_t::face2Minus:   // bottom
      rankNoEnd = size[0] * size[1];
      break;
    case Mesh::face_t::face2Plus:   // top
      rankNoBegin = size[0] * size[1] * (size[2]-1);
      break;
    default:
      break;
    }

    VLOG(1) << "nNumbersTotal: " << nNumbersTotal << ", nDofsPerNode: " << nDofsPerNode
      << "rankNoBegin: " << rankNoBegin << ", rankNoEnd: " << rankNoEnd << ", rankNoStride: " << rankNoStride;

    int i = 0;
    for (int rankNo = rankNoBegin; rankNo < rankNoEnd; rankNo += rankNoStride)
    {
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++, i++)
      {
        assert(i < nNumbersTotal*nDofsPerNode);
        surfaceNumbers[i] = rankNo * nDofsPerNode + nodalDofIndex;
      }
    }
  }
}

}  // namespace
