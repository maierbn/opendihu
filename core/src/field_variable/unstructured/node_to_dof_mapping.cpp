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
  //assert (nodeGlobalNo < (int)nodeDofInformation_.size());
  return nodeDofInformation_[nodeGlobalNo];
}

std::vector<dof_no_t> &NodeToDofMapping::getNodeDofs(node_no_t nodeGlobalNo)
{
  assert (nodeGlobalNo < (int)nodeDofInformation_.size());
  return nodeDofInformation_[nodeGlobalNo].dofs;
}

std::vector<double> &NodeToDofMapping::getNodeScaleFactors(node_no_t nodeGlobalNo)
{
  assert (nodeGlobalNo < (int)nodeDofInformation_.size());
  return nodeDofInformation_[nodeGlobalNo].scaleFactors;
}

node_no_t NodeToDofMapping::nNodes() const
{
  return nodeDofInformation_.size();
}

int NodeToDofMapping::nVersions(node_no_t nodeGlobalNo)
{
  assert (nodeGlobalNo < (int)nodeDofInformation_.size());
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