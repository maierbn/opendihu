#include "field_variable/unstructured/node_to_dof_mapping.h"

#include <cassert>

#include "utility/string_utility.h"
#include "utility/math_utility.h"
#include "easylogging++.h"

namespace FieldVariable
{

bool NodeToDofMapping::containsNode(node_no_t nodeGlobalNo) const
{
  return nodeDofInformation_.find(nodeGlobalNo) != nodeDofInformation_.end();
}

NodeToDofMapping::NodeDofInformation& NodeToDofMapping::getNodeDofInformation(node_no_t nodeGlobalNo)
{
  // this potentially creates the nodeToDof information
  assert (nodeGlobalNo >= 0);
  return nodeDofInformation_[nodeGlobalNo];
}

std::vector<dof_no_t> &NodeToDofMapping::getNodeDofs(node_no_t nodeGlobalNo)
{
  assert (nodeDofInformation_.find(nodeGlobalNo) != nodeDofInformation_.end());
  return nodeDofInformation_[nodeGlobalNo].dofs;
}

dof_no_t NodeToDofMapping::getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  assert (nodeDofInformation_.find(nodeGlobalNo) != nodeDofInformation_.end());
  assert (dofIndex >= 0);
  assert (dofIndex < nodeDofInformation_.at(nodeGlobalNo).dofs.size());
  return nodeDofInformation_.at(nodeGlobalNo).dofs[dofIndex];
}

std::vector<double> &NodeToDofMapping::getNodeScaleFactors(node_no_t nodeGlobalNo)
{
  assert (nodeDofInformation_.find(nodeGlobalNo) != nodeDofInformation_.end());
  return nodeDofInformation_[nodeGlobalNo].scaleFactors;
}

node_no_t NodeToDofMapping::nLocalNodes() const
{
  return nodeDofInformation_.size();
}

int NodeToDofMapping::nVersions(node_no_t nodeGlobalNo)
{
  VLOG(1) << "  NodeToDofMapping::nVersions node " << nodeGlobalNo << ", max: " << nodeDofInformation_.size();
  
  assert (nodeDofInformation_.find(nodeGlobalNo) != nodeDofInformation_.end());
  return nodeDofInformation_[nodeGlobalNo].elementsOfVersion.size();
}

void NodeToDofMapping::output(std::ostream &stream) const
{
  for (auto nodeDofInformation : nodeDofInformation_)
  {
    stream << "[" << nodeDofInformation.first << ":";
    for (int i=0; i<nodeDofInformation.second.dofs.size(); i++)
      stream << nodeDofInformation.second.dofs[i] << " ";
    stream << "] ";
  }
}

std::ostream &operator<<(std::ostream &stream, const NodeToDofMapping &rhs)
{
  rhs.output(stream);
  return stream;
}

std::ostream &operator<<(std::ostream &stream, const NodeToDofMapping::NodeDofInformation::ElementLocalNode &rhs)
{
  stream << rhs.elementGlobalNo << "." << rhs.nodeIdx;
  return stream;
}

};