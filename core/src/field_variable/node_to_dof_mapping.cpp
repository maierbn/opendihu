#include "field_variable/node_to_dof_mapping.h"

#include <cassert>

#include "utility/string_utility.h"
#include "utility/math_utility.h"
#include "easylogging++.h"

namespace FieldVariable
{

bool NodeToDofMapping::containsNode(node_idx_t nodeGlobalNo) const
{
  return nodeDofInformation_.find(nodeGlobalNo) != nodeDofInformation_.end();
}

NodeToDofMapping::NodeDofInformation& NodeToDofMapping::getNodeDofInformation(node_idx_t nodeGlobalNo)
{
  //assert (nodeGlobalNo < (int)nodeDofInformation_.size());
  return nodeDofInformation_[nodeGlobalNo];
}

std::vector<int> &NodeToDofMapping::getNodeDofs(node_idx_t nodeGlobalNo)
{
  assert (nodeGlobalNo < (int)nodeDofInformation_.size());
  return nodeDofInformation_[nodeGlobalNo].dofs;
}

std::vector<double> &NodeToDofMapping::getNodeScaleFactors(node_idx_t nodeGlobalNo)
{
  assert (nodeGlobalNo < (int)nodeDofInformation_.size());
  return nodeDofInformation_[nodeGlobalNo].scaleFactors;
}

node_idx_t NodeToDofMapping::nNodes() const
{
  return nodeDofInformation_.size();
}

int NodeToDofMapping::nVersions(node_idx_t nodeGlobalNo)
{
  assert (nodeGlobalNo < (int)nodeDofInformation_.size());
  return nodeDofInformation_[nodeGlobalNo].elementsOfVersion.size();
}

std::ostream &operator<<(std::ostream &stream, const NodeToDofMapping::NodeDofInformation::ElementLocalNode &rhs)
{
  stream << rhs.elementGlobalNo << "." << rhs.nodeIdx;
  return stream;
}

};