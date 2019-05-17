#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
createNeumannBoundaryConditions(std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<FunctionSpaceType,Quadrature::Gauss<3>, 1>> &neumannBoundaryConditions)
{
  typedef SpatialDiscretization::NeumannBoundaryConditions<FunctionSpaceType,Quadrature::Gauss<3>,1> NeumannBoundaryConditionsType;

  std::array<int,3> nElementsPerCoordinateDirectionLocal;
  for (int i = 0; i < 3; i++)
  {
    nElementsPerCoordinateDirectionLocal[i] = this->meshPartition_->nElementsLocal(i);
  }

  // create neumann BC object
  neumannBoundaryConditions = std::make_shared<NeumannBoundaryConditionsType>(this->context_);

  // determine surface of top and bottom faces of the muscle volume
  double surfaceTopLocal = 0;
  double surfaceBottomLocal = 0;

  // only if the current domain is at the bottom of the overall domain
  if (meshPartition_->ownRankPartitioningIndex(2) == 0)
  {
    // loop over bottom elements
    for (element_no_t elementIndexX = 0; elementIndexX < nElementsPerCoordinateDirectionLocal[0]; elementIndexX++)
    {
      for (element_no_t elementIndexY = 0; elementIndexY < nElementsPerCoordinateDirectionLocal[1]; elementIndexY++)
      {
        element_no_t elementNoLocal = elementIndexY*nElementsPerCoordinateDirectionLocal[1] + elementIndexX;

        // get geometry values
        std::array<Vec3,FunctionSpaceType::nDofsPerElement()> geometryValues;
        this->functionSpace_->getElementGeometry(elementNoLocal, geometryValues);

        if (std::is_same<BasisFunctionType,BasisFunction::LagrangeOfOrder<1>>::value)
        {
          // 2 3
          // 0 1
          Vec3 a = -geometryValues[0] + geometryValues[3];
          Vec3 b = -geometryValues[1] + geometryValues[2];

          // compute surface
          surfaceBottomLocal += 0.5*MathUtility::norm<3>(MathUtility::cross(a,b));
        }
        else if (std::is_same<BasisFunctionType,BasisFunction::LagrangeOfOrder<2>>::value)
        {
          // 6 7 8
          // 3 4 5
          // 0 1 2
          Vec3 a = -geometryValues[0] + geometryValues[8];
          Vec3 b = -geometryValues[2] + geometryValues[6];

          // compute surface
          surfaceBottomLocal += 0.5*MathUtility::norm<3>(MathUtility::cross(a,b));
        }
      }
    }
  }

  // only if the current domain is at the top of the overall domain
  if (meshPartition_->ownRankPartitioningIndex(2) == meshPartition_->nRanks(2)-1)
  {

    // loop over top elements
    element_no_t elementIndexZ = nElementsPerCoordinateDirectionLocal[2]-1;
    for (element_no_t elementIndexX = 0; elementIndexX < nElementsPerCoordinateDirectionLocal[0]; elementIndexX++)
    {
      for (element_no_t elementIndexY = 0; elementIndexY < nElementsPerCoordinateDirectionLocal[1]; elementIndexY++)
      {
        element_no_t elementNoLocal = elementIndexZ*nElementsPerCoordinateDirectionLocal[0]*nElementsPerCoordinateDirectionLocal[1]
          + elementIndexY*nElementsPerCoordinateDirectionLocal[1] + elementIndexX;

        // get geometry values
        std::array<Vec3,FunctionSpaceType::nDofsPerElement()> geometryValues;
        this->functionSpace_->getElementGeometry(elementNoLocal, geometryValues);

        if (std::is_same<BasisFunctionType,BasisFunction::LagrangeOfOrder<1>>::value)
        {
          // 2 3  6 7
          // 0 1  4 5
          Vec3 a = -geometryValues[4] + geometryValues[7];
          Vec3 b = -geometryValues[5] + geometryValues[6];

          // compute surface
          surfaceTopLocal += 0.5*MathUtility::norm<3>(MathUtility::cross(a,b));
        }
        else if (std::is_same<BasisFunctionType,BasisFunction::LagrangeOfOrder<2>>::value)
        {
          // 6 7 8  15 16 17    24 25 26
          // 3 4 5  12 13 14    21 22 23
          // 0 1 2   9 10 11    18 19 20
          Vec3 a = -geometryValues[18] + geometryValues[26];
          Vec3 b = -geometryValues[20] + geometryValues[24];

          // compute surface
          surfaceTopLocal += 0.5*MathUtility::norm<3>(MathUtility::cross(a,b));
          LOG(DEBUG) << "surfaceTopLocal += 0.5*|" << a << " x " << b << "| = " << 0.5*MathUtility::norm<3>(MathUtility::cross(a,b)) << " g: " << geometryValues;
        }
      }
    }
  }

  // reduce global surface values
  double surfaceTopGlobal = 0;
  double surfaceBottomGlobal = 0;

  MPIUtility::handleReturnValue(MPI_Allreduce(&surfaceBottomLocal, &surfaceBottomGlobal, 1, MPI_DOUBLE,
                                              MPI_SUM, this->currentRankSubset_->mpiCommunicator()));
  MPIUtility::handleReturnValue(MPI_Allreduce(&surfaceTopLocal, &surfaceTopGlobal, 1, MPI_DOUBLE,
                                              MPI_SUM, this->currentRankSubset_->mpiCommunicator()));

  // compute flux values that will be set as Neumann BC values
  double fluxTop = 1./surfaceTopGlobal;
  double fluxBottom = -1./surfaceBottomGlobal;

  LOG(DEBUG) << "surfaceTop local: " << surfaceTopLocal << ", global: " << surfaceTopGlobal << ", fluxTop: " << fluxTop;
  LOG(DEBUG) << "surfaceBottom local: " << surfaceBottomLocal << ", global: " << surfaceBottomGlobal << ", fluxBottom: " << fluxBottom;

  // create neumann boundary conditions for elements
  /*
  struct ElementWithFaces
  {
    element_no_t elementNoLocal;   //< the local no of the element

    Mesh::face_t face;              //< face on which the Neumann BC is applied
    std::vector<std::pair<dof_no_t, VecD<nComponents>>> dofVectors;  //< <surface-local dof no, value>, nComponents == FunctionSpaceType::dim() for traction boundary condition or nComponents = 1 for flux BC
    std::vector<dof_no_t> surfaceDofs;    //< dof nos of the volume element that correspond to the face / surface. These are different from the dofs in dofsVector which are numbered for the surface only, surfaceDofs are in the numbering of the volume element.
    // note, for flux BC, dofVectors[i].second is a VecD<1>
  };*/

  typedef typename NeumannBoundaryConditionsType::ElementWithFaces ElementWithFaces;

  std::vector<ElementWithFaces> elementsWithFaces;

  // only if the current domain is at the bottom of the overall domain
  if (meshPartition_->ownRankPartitioningIndex(2) == 0)
  {
    // loop over bottom elements
    for (element_no_t elementIndexX = 0; elementIndexX < nElementsPerCoordinateDirectionLocal[0]; elementIndexX++)
    {
      for (element_no_t elementIndexY = 0; elementIndexY < nElementsPerCoordinateDirectionLocal[1]; elementIndexY++)
      {
        element_no_t elementNoLocal = elementIndexY*nElementsPerCoordinateDirectionLocal[1] + elementIndexX;

        ElementWithFaces elementWithFaces;
        elementWithFaces.elementNoLocal = elementNoLocal;
        elementWithFaces.face = Mesh::face_t::face2Minus;   // facing to bottom

        // add values
        for (int i = 0; i < 4; i++)
        {
          elementWithFaces.dofVectors.push_back(std::pair<dof_no_t, VecD<1>>(i,VecD<1>({fluxBottom})));
        }

        // add volume surface dofs
        // get dofs indices within the numbering of the volume element that correspond to the selected face
        const int D = FunctionSpaceType::dim();
        const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
        const int nVolumeDofsBorder = FunctionSpace::FunctionSpaceBaseDim<D-1,typename FunctionSpaceType::BasisFunction>::nNodesPerElement() * nDofsPerNode;
        std::array<dof_no_t,nVolumeDofsBorder> dofIndices;
        FunctionSpaceType::getFaceDofs(elementWithFaces.face, dofIndices);

        for (int i = 0; i < nVolumeDofsBorder; i++)
        {
          elementWithFaces.surfaceDofs.push_back(dofIndices[i]);
        }

        elementsWithFaces.push_back(elementWithFaces);
      }
    }
  }

  // only if the current domain is at the top of the overall domain
  if (meshPartition_->ownRankPartitioningIndex(2) == meshPartition_->nRanks(2)-1)
  {
    element_no_t elementIndexZ = nElementsPerCoordinateDirectionLocal[2]-1;
    for (element_no_t elementIndexX = 0; elementIndexX < nElementsPerCoordinateDirectionLocal[0]; elementIndexX++)
    {
      for (element_no_t elementIndexY = 0; elementIndexY < nElementsPerCoordinateDirectionLocal[1]; elementIndexY++)
      {
        element_no_t elementNoLocal = elementIndexZ*nElementsPerCoordinateDirectionLocal[0]*nElementsPerCoordinateDirectionLocal[1]
          + elementIndexY*nElementsPerCoordinateDirectionLocal[1] + elementIndexX;

        ElementWithFaces elementWithFaces;
        elementWithFaces.elementNoLocal = elementNoLocal;
        elementWithFaces.face = Mesh::face_t::face2Plus;   // facing to top

        // add values
        for (int i = 0; i < 4; i++)
        {
          elementWithFaces.dofVectors.push_back(std::pair<dof_no_t, VecD<1>>(i,VecD<1>({fluxTop})));
        }

        // add volume surface dofs
        // get dofs indices within the numbering of the volume element that correspond to the selected face
        const int D = FunctionSpaceType::dim();
        const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
        const int nVolumeDofsBorder = FunctionSpace::FunctionSpaceBaseDim<D-1,typename FunctionSpaceType::BasisFunction>::nNodesPerElement() * nDofsPerNode;
        std::array<dof_no_t,nVolumeDofsBorder> dofIndices;
        FunctionSpaceType::getFaceDofs(elementWithFaces.face, dofIndices);

        for (int i = 0; i < nVolumeDofsBorder; i++)
        {
          elementWithFaces.surfaceDofs.push_back(dofIndices[i]);
        }

        elementsWithFaces.push_back(elementWithFaces);
      }
    }
  }

  neumannBoundaryConditions->initialize(this->functionSpace_, elementsWithFaces);
}

} // namespace
