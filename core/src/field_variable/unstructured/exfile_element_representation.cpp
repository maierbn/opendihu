#include "field_variable/unstructured/exfile_element_representation.h"

#include <sstream>
#include <fstream>
#include <cassert>

#include "utility/string_utility.h"
#include "easylogging++.h"
#include "utility/math_utility.h"

namespace FieldVariable
{

using namespace StringUtility;
  
void ExfileElementRepresentation::parseFromExelemFile(std::string content)
{
  VLOG(1) << "ExfileElementRepresentation::parseFromExelemFile(" << content << ")";
  
  int nNodes = 0;
  int nodeNo = 0;
  int nValues = 0;
  std::stringstream nodeNoStr;
  
  size_t pos = 0;
  while(pos < content.size())
  {
    // extract next line
    size_t posNewline = content.find("\n",pos);
    std::string line = content.substr(pos, posNewline-pos);
    if (posNewline == std::string::npos)
      pos = content.size();
    else
      pos = posNewline+1;
   
    
    //VLOG(1) << " line [" << line << "]";
    
    if (line.find("#Nodes=") != std::string::npos)
    {
      nNodes = getNumberAfterString(line, "#Nodes=");
      node_.resize(nNodes);
      nodeNo = 0;
      nodeNoStr.str("");
      nodeNoStr << nodeNo+1 << ".";
      VLOG(1) << "nNodes: " << nNodes;
    }
    else
    {
      if (line.find(nodeNoStr.str()) != std::string::npos)
      {
        nValues = getNumberAfterString(line, "#Values=");
          VLOG(1) << " nValues: " << nValues;
        node_[nodeNo].valueIndices.resize(nValues);
        node_[nodeNo].scaleFactorIndices.resize(nValues);
        
        // prepare next block
        nodeNo++;
        nodeNoStr.str("");
        nodeNoStr << nodeNo+1 << ".";
      }
      else if (line.find("Value indices:") != std::string::npos)
      {
        extractUntil(line, "Value indices:");
        for(int i=0; i<nValues; i++)
        {
          trim(line);
          node_[nodeNo-1].valueIndices[i] = atoi(line.c_str())-1;
          extractUntil(line, " ");
          //VLOG(1) << " i=" << i << ", valueIndices: " << node_[nodeNo].valueIndices[i];
        }
      }
      else if (line.find("Scale factor indices:") != std::string::npos)
      {
        extractUntil(line, "Scale factor indices:");
        for(int i=0; i<nValues; i++)
        {
          trim(line);
          node_[nodeNo-1].scaleFactorIndices[i] = atoi(line.c_str())-1;
          extractUntil(line, " ");
          //VLOG(1) << " i=" << i << ", scaleFactorIndices: " << node_[nodeNo].scaleFactorIndices[i];
        }
      }
    }
  }
}

void ExfileElementRepresentation::setNumberNodes(int nNodes)
{
  node_.resize(nNodes);
}

bool ExfileElementRepresentation::operator==(const ExfileElementRepresentation& rhs) const
{
  //VLOG(1) << "    exfileElementRepresentation sizes: " << node_.size() << ", " << rhs.node_.size();
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

void ExfileElementRepresentation::outputHeaderExelemFile(std::ostream &file)
{
  int no = 1;
  for (std::vector<Node>::const_iterator iter = node_.cbegin(); iter != node_.cend(); iter++, no++)
  {
    file << "   " << no << ". #Values=" << iter->valueIndices.size() << std::endl
      << "     Value indices:";
    
    for (unsigned int i=0; i<iter->valueIndices.size(); i++)
    {
      file << " " << iter->valueIndices[i]+1;
    }
    file << std::endl
      << "     Scale factor indices:";
    
    unsigned int i=0;
    for (; i<iter->scaleFactorIndices.size(); i++)
    {
      file << " " << iter->scaleFactorIndices[i]+1;
    }
    // fill the rest of the scale factors with 0s (means scale factor of 1.0)
    for (; i<iter->valueIndices.size(); i++)
    {
      file << " 0";
    }
    file << std::endl;
  }
}

void ExfileElementRepresentation::output(std::ostream &stream) const
{
  stream << "ExfileElementRepresentation: " << node_.size() << " nodes: ";
  for (unsigned int i=0; i<node_.size(); i++)
  {
    stream << i << ":(";
    for (auto valueIndex : node_[i].valueIndices)
      stream << valueIndex+1 << ",";
    stream << ") ";
  }
}

std::ostream &operator<<(std::ostream &stream, const ExfileElementRepresentation &rhs)
{
  rhs.output(stream);
  return stream; 
}

};
