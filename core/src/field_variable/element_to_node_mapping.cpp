#include "field_variable/element_to_node_mapping.h"

#include <cassert>

#include "utility/string_utility.h"
#include "easylogging++.h"

namespace FieldVariable
{
  
using namespace StringUtility;
  
void ElementToNodeMapping::setNumberElements(int nElements)
{
  elements_.resize(nElements);
}

void ElementToNodeMapping::parseElementFromExelemFile(std::string content)
{
  int elementNo = 0;
  bool nodesFollow = false;
  bool scaleFactorsFollow = false;
  
  int pos = 0;
  while(pos < (int)content.size())
  {
    // extract next line
    unsigned int posNewline = content.find("\n",pos);
    std::string line = content.substr(pos, posNewline);
    if (posNewline == std::string::npos)
      pos = content.size();
    else
      pos = posNewline+1;
    
    if (line.find("Element:") != std::string::npos)
    {
      elementNo = getNumberAfterString(line,"Element:");
    }
    
    if (line.find("Nodes:") != std::string::npos)
    {
      nodesFollow = true;
      continue;
    }
    
    if (line.find(" Scale factors:") != std::string::npos)
    {
      scaleFactorsFollow = true;
      continue;
    }
    
    if (nodesFollow)
    {
      trim(line);
      while(!line.empty())
      {
        int globalNodeNo = atoi(line.c_str())-1;
        elements_[elementNo].globalNodeNo.push_back(globalNodeNo);
        extractUntil(line, " ");
        trim(line);
      }
      nodesFollow = false;
    }
    
    if (scaleFactorsFollow)
    {
      trim(line);
      while(!line.empty())
      {
        double scaleFactor = atof(line.c_str());
        elements_[elementNo].scaleFactors.push_back(scaleFactor);
        extractUntil(line, " ");
        trim(line);
      }
    }
  }
}

ElementToNodeMapping::Element& ElementToNodeMapping::getElement(int elementNo)
{
  if (elementNo >= (int)elements_.size())
  {
    LOG(ERROR) << "element no out of range, " << elementNo << ">=" << elements_.size();
  }
  return elements_[elementNo];
}

void ElementToNodeMapping::outputElementExelemFile(std::ofstream &file, element_idx_t elementGlobalNo)
{
  assert(elementGlobalNo < (int)elements_.size());
  
  file << " Element: " << elementGlobalNo+1 << " 0 0" << std::endl 
    << " Nodes:" << std::endl;
    
  // output global node numbers
  for (unsigned int nodeNo = 0; nodeNo < elements_[elementGlobalNo].globalNodeNo.size(); nodeNo++)
  {
    file << " " << elements_[elementGlobalNo].globalNodeNo[nodeNo]+1;
  }
  file << std::endl 
    << "Scale factors:" << std::endl << " ";
    
  // output scale factors
  StringUtility::outputValuesBlock(file, elements_[elementGlobalNo].scaleFactors, 5);
}
  
};