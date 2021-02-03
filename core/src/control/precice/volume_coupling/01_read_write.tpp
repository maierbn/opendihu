#include "control/precice/volume_coupling/01_read_write.h"

#include <sstream>

#include "spatial_discretization/neumann_boundary_conditions/00_neumann_boundary_conditions_base.h"
#include "utility/vector_operators.h"

namespace Control
{

#ifdef HAVE_PRECICE
template<typename NestedSolver>
void PreciceAdapterVolumeCouplingReadWrite<NestedSolver>::
preciceReadData()
{
  if (!this->preciceSolverInterface_->isReadDataAvailable())
    return;

  LOG(DEBUG) << "read data from precice";

  using SlotConnectorDataType = typename NestedSolver::Data::SlotConnectorDataType;
  std::shared_ptr<SlotConnectorDataType> slotConnectorData = this->nestedSolver_.data().getSlotConnectorData();

  int nArrayItems = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(slotConnectorData);   // number of fibers if there are fibers

  // loop over data
  for (typename PreciceAdapterVolumeCouplingInitialize<NestedSolver>::PreciceData &preciceData : this->preciceData_)
  {
    if (preciceData.ioType == PreciceAdapterVolumeCouplingReadWrite<NestedSolver>::PreciceData::ioRead)
    {
      int nEntries = preciceData.preciceMesh->nNodesLocal;

      if (preciceData.isGeometryField)
      {
        nEntries = preciceData.preciceMesh->nNodesLocal*3;
      }

      // allocate temporary memory
      scalarValues_.resize(nEntries);

      // get all data at once
      this->preciceSolverInterface_->readBlockVectorData(preciceData.preciceDataId, preciceData.preciceMesh->nNodesLocal,
                                                         preciceData.preciceMesh->preciceVertexIds.data(), scalarValues_.data());

      // get the mesh partition
      std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase
        = SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase(slotConnectorData, preciceData.slotNo, 0);

      int nDofsLocalWithoutGhosts = meshPartitionBase->nDofsLocalWithoutGhosts();

      // get the vector of values [0,1,...,nDofsLocalWithGhosts]
      const std::vector<PetscInt> &dofNosLocalWithGhosts = meshPartitionBase->dofNosLocal();
      std::vector<PetscInt> dofNosLocalWithoutGhosts(dofNosLocalWithGhosts.begin(), dofNosLocalWithGhosts.begin()+nDofsLocalWithoutGhosts);

      // store received data in field variable
      if (preciceData.isGeometryField)
      {
        // loop over fibers if there are any
        for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++)
        {
          // fill the vector geometryValues_ with the geometry values of the current fiber or mesh
          geometryValues_.resize(nDofsLocalWithoutGhosts);
          for (int dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts; dofNoLocal++)
          {
            for (int componentNo = 0; componentNo < 3; componentNo++)
            {
              geometryValues_[dofNoLocal][componentNo] = scalarValues_[3*(arrayIndex*nDofsLocalWithoutGhosts + dofNoLocal) + componentNo];
            }
          }

          SlotConnectorDataHelper<SlotConnectorDataType>::slotSetGeometryValues(slotConnectorData, preciceData.slotNo, arrayIndex,
                                                                                dofNosLocalWithoutGhosts, geometryValues_);
        }
      }
      else
      {
        // loop over fibers if there are any
        for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++)
        {
          // fill the vector geometryValues_ with the geometry values of the current fiber or mesh
          scalarValuesOfMesh_.assign(scalarValues_.begin() +  arrayIndex   *nDofsLocalWithoutGhosts,
                                     scalarValues_.begin() + (arrayIndex+1)*nDofsLocalWithoutGhosts);

          SlotConnectorDataHelper<SlotConnectorDataType>::slotSetValues(slotConnectorData, preciceData.slotNo, arrayIndex,
                                                                        dofNosLocalWithoutGhosts, scalarValuesOfMesh_);
        }
      }
    }
  }
}

template<typename NestedSolver>
void PreciceAdapterVolumeCouplingReadWrite<NestedSolver>::
preciceWriteData()
{
  if (!this->preciceSolverInterface_->isWriteDataRequired(this->timeStepWidth_))
    return;

  // write data to precice
  LOG(DEBUG) << "write data to precice";

  using SlotConnectorDataType = typename NestedSolver::Data::SlotConnectorDataType;
  std::shared_ptr<SlotConnectorDataType> slotConnectorData = this->nestedSolver_.data().getSlotConnectorData();

  int nArrayItems = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(slotConnectorData);   // number of fibers if there are fibers

  // loop over data
  for (typename PreciceAdapterVolumeCouplingInitialize<NestedSolver>::PreciceData &preciceData : this->preciceData_)
  {
    if (preciceData.ioType == PreciceAdapterVolumeCouplingInitialize<NestedSolver>::PreciceData::ioWrite)
    {
      scalarValues_.clear();

      // get the mesh partition
      std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase
        = SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase(slotConnectorData, preciceData.slotNo, 0);

      int nDofsLocalWithoutGhosts = meshPartitionBase->nDofsLocalWithoutGhosts();

      // get the vector of values [0,1,...,nDofsLocalWithGhosts]
      const std::vector<PetscInt> &dofNosLocalWithGhosts = meshPartitionBase->dofNosLocal();
      std::vector<PetscInt> dofNosLocalWithoutGhosts(dofNosLocalWithGhosts.begin(), dofNosLocalWithGhosts.begin()+nDofsLocalWithoutGhosts);

      // if it is a geometry field, get the node positions of a mesh
      if (preciceData.isGeometryField)
      {
        geometryValues_.clear();

        // loop over fibers if there are any
        for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++)
        {
          SlotConnectorDataHelper<SlotConnectorDataType>::slotGetGeometryValues(slotConnectorData, preciceData.slotNo, arrayIndex,
                                                                                dofNosLocalWithoutGhosts, geometryValues_);
        }
        LOG(DEBUG) << "Using geometry field of opendihu meshes of " << nArrayItems << " array items (e.g., fibers).";

        // transform to contiguous memory layout for precice
        scalarValues_.resize(3*geometryValues_.size());

        for (int entryNo = 0; entryNo < geometryValues_.size(); entryNo++)
          for (int componentNo = 0; componentNo < 3; componentNo++)
            scalarValues_[3*entryNo + componentNo] = geometryValues_[entryNo][componentNo];

      }
      else
      {
        // the data is a normal slot, no geometry field

        // loop over fibers if there are any
        for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++)
        {
          // static void slotGetValues(std::shared_ptr<SlotConnectorDataType> slotConnectorData,
          //   int slotNo, int arrayIndex, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values);
          SlotConnectorDataHelper<SlotConnectorDataType>::slotGetValues(slotConnectorData, preciceData.slotNo, arrayIndex,
                                                                        dofNosLocalWithoutGhosts, scalarValues_);
        }
      }

#ifndef NDEBUG
      LOG(DEBUG) << "write data to precice: " << scalarValues_;
#endif
      // scale the values by a factor given in the config
      for (double &value : scalarValues_)
        value *= this->scalingFactor_;

      // write values in precice
      this->preciceSolverInterface_->writeBlockVectorData(preciceData.preciceDataIdDisplacements, preciceData.preciceMesh->nNodesLocal,
                                                          preciceData.preciceMesh->preciceVertexIds.data(), scalarValues_.data());
    }
  }

  LOG(DEBUG) << "write data to precice complete";
}
#endif

}  // namespace
