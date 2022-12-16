#include "control/precice/surface_coupling/02_read_write.h"

#include <sstream>

#include "spatial_discretization/neumann_boundary_conditions/00_neumann_boundary_conditions_base.h"
#include "utility/vector_operators.h"

namespace Control
{

#ifdef HAVE_PRECICE
template<typename NestedSolver>
void PreciceAdapterReadWrite<NestedSolver>::
preciceReadData()
{
  if (!this->preciceSolverInterface_->isReadDataAvailable())
    return;

  LOG(DEBUG) << "read data from precice";

  // loop over data
  for (typename PreciceAdapterInitialize<NestedSolver>::PreciceData &preciceData : this->preciceData_)
  {
    if (preciceData.ioType == PreciceAdapterReadWrite<NestedSolver>::PreciceData::ioRead)
    {
      // allocate memory
      int nEntries = preciceData.preciceMesh->nNodesLocal*3;

      // if the data is displacements and velocities
      if (!preciceData.displacementsName.empty())
      {
        displacementValues_.resize(nEntries);
        velocityValues_.resize(nEntries);

        // get all data at once
        this->preciceSolverInterface_->readBlockVectorData(preciceData.preciceDataIdDisplacements, preciceData.preciceMesh->nNodesLocal,
                                                           preciceData.preciceMesh->preciceVertexIds.data(), displacementValues_.data());

        // this->preciceSolverInterface_->readBlockVectorData(preciceData.preciceDataIdVelocities, preciceData.preciceMesh->nNodesLocal,
        //                                                    preciceData.preciceMesh->preciceVertexIds.data(), velocityValues_.data());

        setDirichletBoundaryConditions(preciceData);
      }
      // if the data is traction
      else if (!preciceData.tractionName.empty())
      {
        tractionValues_.resize(nEntries);
        this->preciceSolverInterface_->readBlockVectorData(preciceData.preciceDataIdTraction, preciceData.preciceMesh->nNodesLocal,
                                                           preciceData.preciceMesh->preciceVertexIds.data(), tractionValues_.data());

        setNeumannBoundaryConditions(preciceData);
      }
      else
      {
        LOG(FATAL) << "Unknown precice data (read), none of displacements, velocities or traction is set.";
      }
    }
  }
}

template<typename NestedSolver>
void PreciceAdapterReadWrite<NestedSolver>::
setDirichletBoundaryConditions(typename PreciceAdapterInitialize<NestedSolver>::PreciceData &preciceData)
{
  std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBCValues;

  // if this rank has no data, do not set any boundary conditions
  if (preciceData.preciceMesh->nNodesLocal != 0)
  {
    // loop over nodes
    const int nNodesX = this->functionSpace_->nNodesLocalWithoutGhosts(0);
    const int nNodesY = this->functionSpace_->nNodesLocalWithoutGhosts(1);
    const int nNodesZ = this->functionSpace_->nNodesLocalWithoutGhosts(2);

    // set node index in z direction for bottom surface
    int nodeIndexZ = 0;

    // for top surface
    if (preciceData.preciceMesh->face == PreciceAdapterReadWrite<NestedSolver>::PreciceMesh::face2Plus)
    {
      nodeIndexZ = nNodesZ-1;
    }

    // loop over nodes to set the received values
    newDirichletBCValues.reserve(nNodesX * nNodesY);
    int valueIndex = 0;

    // loop over nodes of surface mesh
    for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++)
    {
      for (int nodeIndexX = 0; nodeIndexX < nNodesX; nodeIndexX++, valueIndex++)
      {
        node_no_t nodeNoLocal =
          nodeIndexZ * nNodesX * nNodesY
          + nodeIndexY * nNodesX
          + nodeIndexX;

        dof_no_t dofNoLocal = nodeNoLocal;
        global_no_t dofNoGlobal = this->functionSpace_->meshPartition()->getDofNoGlobalPetsc(dofNoLocal);

        // assign received values to dirichlet bc vector of size 6
        std::array<double,6> newDirichletBCValue;

        for (int i = 0; i < 3; i++)
        {
          newDirichletBCValue[i] = displacementValues_[3*valueIndex + i];
          newDirichletBCValue[3+i] = velocityValues_[3*valueIndex + i];
        }

        newDirichletBCValues.push_back(std::pair<global_no_t,std::array<double,6>>(dofNoGlobal, newDirichletBCValue));
      }
    }

    LOG(INFO) << "read data from precice complete, displacement values: " << displacementValues_
      << ", velocityValues: " << velocityValues_;
    LOG(INFO) << "read and set Dirichlet BC: " << newDirichletBCValues;
  }

  //! set new dirichlet boundary condition values
  this->updateDirichletBoundaryConditions(this->nestedSolver_, newDirichletBCValues);
}

template<typename NestedSolver>
void PreciceAdapterReadWrite<NestedSolver>::
setNeumannBoundaryConditions(typename PreciceAdapterInitialize<NestedSolver>::PreciceData &preciceData)
{
  // set traction values as neumann boundary conditions
  using FunctionSpace = typename PreciceAdapterNestedSolver<NestedSolver>::FunctionSpace;
  using ElementWithFacesType = typename SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,Quadrature::Gauss<3>,3>::ElementWithFaces;
  /*
  struct ElementWithFaces
  {
    element_no_t elementNoLocal;                                     //< the local no of the element

    Mesh::face_t face;                                               //< face on which the Neumann BC is applied
    std::vector<std::pair<dof_no_t, VecD<nComponents>>> dofVectors;  //< <surface-local dof no, value>, nComponents == FunctionSpaceType::dim() for traction boundary condition or nComponents = 1 for flux BC
    std::vector<dof_no_t> surfaceDofs;                               //< dof nos of the volume element that correspond to the face / surface. These are different from the dofs in dofsVector which are numbered for the surface only, surfaceDofs are in the numbering of the volume element.
    // note, for flux BC, dofVectors[i].second is a VecD<1>
  };*/

  std::vector<ElementWithFacesType> neumannBoundaryConditionElements;

  // if this rank has no data, do not set any boundary conditions
  if (preciceData.preciceMesh->nNodesLocal != 0)
  {

    const int nNodesX = this->functionSpace_->nNodesLocalWithoutGhosts(0);
    const int nNodesY = this->functionSpace_->nNodesLocalWithoutGhosts(1);
    const int nNodesZ = this->functionSpace_->nNodesLocalWithoutGhosts(2);

    int nElementsX = this->functionSpace_->meshPartition()->nElementsLocal(0);
    int nElementsY = this->functionSpace_->meshPartition()->nElementsLocal(1);
    int nElementsZ = this->functionSpace_->meshPartition()->nElementsLocal(2);

    // set node and element indics in z direction for bottom surface
    int nodeIndexZ = 0;
    int elementalNodeIndexZ = 0;
    int elementIndexZ = 0;

    // for top surface
    if (preciceData.preciceMesh->face == PreciceAdapterReadWrite<NestedSolver>::PreciceMesh::face2Plus)
    {
      nodeIndexZ = nNodesZ-1;
      elementIndexZ = nElementsZ-1;
      elementalNodeIndexZ = 2;
    }

    // loop over elements
    for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++)
    {
      for (int elementIndexX = 0; elementIndexX < nElementsX; elementIndexX++)
      {
        ElementWithFacesType elementWithFaces;
        element_no_t elementNoLocal = elementIndexZ * nElementsX*nElementsY + elementIndexY * nElementsX + elementIndexX;
        elementWithFaces.elementNoLocal = elementNoLocal;

        // set surface dofs
        Mesh::face_t face = Mesh::face_t::face2Minus;
        if (preciceData.preciceMesh->face == PreciceAdapterReadWrite<NestedSolver>::PreciceMesh::face2Plus)
          face = Mesh::face_t::face2Plus;

        elementWithFaces.face = face;

        // get dofs indices within the numbering of the volume element that correspond to the selected face
        const int nDofsPerNode = FunctionSpace::nDofsPerNode();
        const int nSurfaceDofs = ::FunctionSpace::FunctionSpaceBaseDim<2,typename FunctionSpace::BasisFunction>::nNodesPerElement() * nDofsPerNode;
        std::array<dof_no_t,nSurfaceDofs> surfaceDofs;
        FunctionSpace::getFaceDofs(face, surfaceDofs);

        elementWithFaces.surfaceDofs.assign(surfaceDofs.begin(), surfaceDofs.end());

        // loop over the nodes of the element
        for (int elementalNodeIndexY = 0; elementalNodeIndexY < 3; elementalNodeIndexY++)
        {
          for (int elementalNodeIndexX = 0; elementalNodeIndexX < 3; elementalNodeIndexX++)
          {
            int elementalDofIndex = elementalNodeIndexZ * 9 + elementalNodeIndexY * 3 + elementalNodeIndexX;

            dof_no_t dofNoLocal = this->functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);

            // only iterate over local dofs, the ghost dofs are not considered here (although it would be correct to communicate them with precice)
            if (dofNoLocal >= this->functionSpace_->nDofsLocalWithoutGhosts())
              continue;

            int valueIndex = dofNoLocal - nodeIndexZ*nNodesX*nNodesY;

            LOG(DEBUG) << "(x,y,z)=(" << elementalNodeIndexX << "," << elementalNodeIndexY << "," << elementalNodeIndexZ << ") dofNoLocal "
              << dofNoLocal << ", valueIndex: " << valueIndex << "/" << tractionValues_.size()/3;

            assert(valueIndex >= 0);
            if (valueIndex >= preciceData.preciceMesh->nNodesLocal)
            {
              LOG(ERROR) << "valueIndex: " << valueIndex << ", dofNoLocal: " << dofNoLocal << ", nNodes: " << nNodesX
                << "," << nNodesY << ", nodeIndexZ: " << nodeIndexZ << ", nNodesLocal: " << preciceData.preciceMesh->nNodesLocal
                << " elementNoLocal: " << elementNoLocal << ", elementalDofIndex: " << elementalDofIndex
                << ", elementalNodeIndexZ: " << elementalNodeIndexZ << ", meshPartition: " << *this->functionSpace_->meshPartition();
            }
            assert(valueIndex < preciceData.preciceMesh->nNodesLocal);


            Vec3 traction;
            for (int i = 0; i < 3; i++)
            {
              traction[i] = tractionValues_[3*valueIndex + i];
            }

            dof_no_t surfaceDof = elementalDofIndex;

            // for top surface
            if (preciceData.preciceMesh->face == PreciceAdapterReadWrite<NestedSolver>::PreciceMesh::face2Plus)
            {
              surfaceDof = 18+elementalDofIndex;
            }
            elementWithFaces.dofVectors.push_back(std::pair<dof_no_t,Vec3>(surfaceDof, traction));

            //LOG(INFO) << "dofVectors: " << elementWithFaces.dofVectors << ", traction: " << traction;
          }
        }
        neumannBoundaryConditionElements.push_back(elementWithFaces);
      }
    }
  }

  // create new Neumann BC object
  using NeumannBoundaryConditionsType = SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,Quadrature::Gauss<3>,3>;
  std::shared_ptr<NeumannBoundaryConditionsType> neumannBoundaryConditions = std::make_shared<NeumannBoundaryConditionsType>(this->context_);
  neumannBoundaryConditions->initialize(this->functionSpace_, neumannBoundaryConditionElements);
  neumannBoundaryConditions->setDeformationGradientField(this->deformationGradientField(this->nestedSolver_));

#ifndef NDEBUG
  std::stringstream s;
  for (int i = 0; i < neumannBoundaryConditionElements.size(); i++)
  {
    s << "{el. " << neumannBoundaryConditionElements[i].elementNoLocal
      << ", \"" << Mesh::getString(neumannBoundaryConditionElements[i].face)
      << "\", dofVectors: [";
    for (int j = 0; j < neumannBoundaryConditionElements[i].dofVectors.size(); j++)
    {
      s << "(" << neumannBoundaryConditionElements[i].dofVectors[j].first
        << ": (";
      for (int k = 0; k < neumannBoundaryConditionElements[i].dofVectors[j].second.size(); k++)
      {
        if (k != 0)
          s << ",";
        s << neumannBoundaryConditionElements[i].dofVectors[j].second[k];
      }
      s << ")),";
    }
    s << "], surfaceDofs: ";
    for (int j = 0; j < neumannBoundaryConditionElements[i].surfaceDofs.size(); j++)
      s << neumannBoundaryConditionElements[i].surfaceDofs[j] << ",";
    s << "} ";
  }
  LOG(DEBUG) << "read and set Neumann BC:\n" << s.str();
#endif

  // set Neumann BCs in the static hyperelasticity of the TimeSteppingScheme::DynamicHyperelasticitySolver solver
  this->updateNeumannBoundaryConditions(this->nestedSolver_, neumannBoundaryConditions);

  LOG(DEBUG) << "read data from precice complete, traction values: " << tractionValues_;
}

template<typename NestedSolver>
void PreciceAdapterReadWrite<NestedSolver>::
preciceWriteData()
{
  if (!this->preciceSolverInterface_->isWriteDataRequired(this->timeStepWidth_))
    return;

  // write data to precice
  LOG(DEBUG) << "write data to precice";

  // loop over data
  for (typename PreciceAdapterInitialize<NestedSolver>::PreciceData &preciceData : this->preciceData_)
  {
    if (preciceData.ioType == PreciceAdapterInitialize<NestedSolver>::PreciceData::ioWrite)
    {
      // if the data is displacements and velocities
      if (!preciceData.displacementsName.empty())
      {
        // convert geometry values to precice data layout
        displacementValues_.clear();
        velocityValues_.clear();

        this->getDisplacementVelocityValues(this->nestedSolver_, preciceData.preciceMesh->dofNosLocal, displacementValues_, velocityValues_);

        LOG(INFO) << "write displacements data to precice: " << displacementValues_;
        //LOG(INFO) << "write velocities data to precice: " << velocityValues_;


        // scale displacement and velocity values
        for (double &value : displacementValues_)
          value *= this->scalingFactor_;

        // for (double &value : velocityValues_)
        //   value *= this->scalingFactor_;

        // write displacement values in precice
        this->preciceSolverInterface_->writeBlockVectorData(preciceData.preciceDataIdDisplacements, preciceData.preciceMesh->nNodesLocal,
                                                            preciceData.preciceMesh->preciceVertexIds.data(), displacementValues_.data());

        // write velocity values in precice
        // this->preciceSolverInterface_->writeBlockVectorData(preciceData.preciceDataIdVelocities, preciceData.preciceMesh->nNodesLocal,
        //                                                     preciceData.preciceMesh->preciceVertexIds.data(), velocityValues_.data());
      }
      // if the data is traction
      else if (!preciceData.tractionName.empty())
      {
        // convert geometry values to precice data layout
        tractionValues_.clear();
        this->getTractionValues(this->nestedSolver_, preciceData.preciceMesh->dofNosLocal, tractionValues_);

#ifndef NDEBUG
        LOG(DEBUG) << "write traction data to precice: " << tractionValues_;
        std::stringstream s;
        for (int i = 2; i < tractionValues_.size(); i+=3)
        {
          s << " " << tractionValues_[i];
        }
        LOG(DEBUG) << "z values of traction: " << s.str();
#endif
        // scale traction values, they are always scaled by the factor of -1
        for (double &value : tractionValues_)
        {
          value *= this->scalingFactor_;
          value *= -1;
        }

        this->preciceSolverInterface_->writeBlockVectorData(preciceData.preciceDataIdTraction, preciceData.preciceMesh->nNodesLocal,
                                                            preciceData.preciceMesh->preciceVertexIds.data(), tractionValues_.data());
      }
      else
      {
        LOG(FATAL) << "Unknown precice data (write), none of displacements, velocities or traction is set.";
      }
    }
  }

  LOG(DEBUG) << "write traction data to precice complete";
}
#endif

}  // namespace
