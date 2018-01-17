#include "field_variable/exfile_element_representation.h"

#include <sstream>
#include <fstream>
#include <cassert>

#include "utility/string_utility.h"
#include "easylogging++.h"

namespace FieldVariable
{

using namespace StringUtility;
  
void ExfileElementRepresentation::parseFromExelemFile(std::string content)
{
  int nNodes = 0;
  int nodeNo = 0;
  int nValues = 0;
  std::stringstream nodeNoStr;
  
  size_t pos = 0;
  while(pos < content.size())
  {
    // extract next line
    size_t posNewline = content.find("\n",pos);
    std::string line = content.substr(pos, posNewline);
    if (posNewline == std::string::npos)
      pos = content.size();
    else
      pos = posNewline+1;
   
    if (line.find("#Nodes=") != std::string::npos)
    {
      nNodes = getNumberAfterString(line, "#Nodes=");
      node_.resize(nNodes);
      nodeNo = 0;
      nodeNoStr.str("");
      nodeNoStr << nodeNo+1 << ".";
    }
    else
    {
      if (line.find(nodeNoStr.str()) != std::string::npos)
      {
        // prepare next block
        nodeNo++;
        nodeNoStr.str("");
        nodeNoStr << nodeNo << ".";
        nValues = getNumberAfterString(line, "#Values=");
        node_[nodeNo].valueIndices.resize(nValues);
        node_[nodeNo].scaleFactorIndices.resize(nValues);
      }
      else if (line.find("Value indices:") != std::string::npos)
      {
        extractUntil(line, "Value indices:");
        for(int i=0; i<nValues; i++)
        {
          trim(line);
          node_[nodeNo].valueIndices[i] = atoi(line.c_str());
          extractUntil(line, " ");
        }
      }
      else if (line.find("Scale factor indices:") != std::string::npos)
      {
        extractUntil(line, "Scale factor indices:");
        for(int i=0; i<nValues; i++)
        {
          trim(line);
          node_[nodeNo].scaleFactorIndices[i] = atoi(line.c_str());
          extractUntil(line, " ");
        }
      }
    }
  }
}

bool ExfileElementRepresentation::operator==(const ExfileElementRepresentation& rhs) const
{
  if (node_.size() != rhs.node_.size())
    return false;
  
  for (size_t i=0; i<node_.size(); i++)
  {
    if (node_[i].valueIndices.size() != rhs.node_[i].valueIndices.size())
      return false;
    if (node_[i].scaleFactorIndices.size() != rhs.node_[i].scaleFactorIndices.size())
      return false;
    
    for (size_t j=0; j<node_[i].valueIndices.size(); j++)
    {
      if (node_[i].valueIndices[j] != rhs.node_[i].valueIndices[j])
        return false;
    }
    for (size_t j=0; j<node_[i].scaleFactorIndices.size(); j++)
    {
      if (node_[i].scaleFactorIndices[j] != rhs.node_[i].scaleFactorIndices[j])
        return false;
    }
  }
  return true;
}

ExfileElementRepresentation::Node& ExfileElementRepresentation::getNode(int nodeIndex)
{
  assert(nodeIndex < (int)node_.size());
  return node_[nodeIndex];
}

void ExfileElementRepresentation::outputHeaderExelemFile(std::ofstream &file)
{
  int no = 1;
  for (std::vector<Node>::const_iterator iter = node_.cbegin(); iter != node_.cend(); iter++, no++)
  {
    file << " " << no << ". #Values=" << iter->valueIndices.size() << std::endl
      << "Value indices:";
    
    for (unsigned int i=0; i<iter->valueIndices.size(); i++)
    {
      file << " " << iter->valueIndices[i];
    }
    file << std::endl
      << "Scale factor indices:";
    
    for (unsigned int i=0; i<iter->valueIndices.size(); i++)
    {
      file << " " << iter->valueIndices[i];
    }
    file << std::endl;
  }
}

};