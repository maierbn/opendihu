#include "field_variable/element_to_node_mapping.h"

#include <cassert>

#include "utility/string_utility.h"
#include "easylogging++.h"

namespace FieldVariable
{
  
using namespace StringUtility;
  
void ElementToNodeMapping::setNumberElements(element_no_t nElements)
{
  elements_.resize(nElements);
}

void ElementToNodeMapping::parseElementFromExelemFile(std::string content)
{
  VLOG(2) << "ElementToNodeMapping::parseElementFromExelemFile(" << content << "), elements_.size: " << elements_.size();
  int elementNo = 0;
  bool nodesFollow = false;
  bool scaleFactorsFollow = false;
  
  int pos = 0;
  while(pos < (int)content.size())
  {
    // extract next line
    unsigned int posNewline = content.find("\n",pos);
    std::string line = content.substr(pos, posNewline-pos);
    if (posNewline == std::string::npos)
      pos = content.size();
    else
      pos = posNewline+1;
    
    //VLOG(2) << "   line [" << StringUtility::replace(line,"\r","") << "],pos: " << pos << ", size: " << content.size();
    
    if (line.find("Element:") != std::string::npos)
    {
      elementNo = getNumberAfterString(line,"Element:")-1;
      VLOG(2) << "   elementNo=" << elementNo;
    }
    
    if (line.find("Nodes:") != std::string::npos)
    {
      nodesFollow = true;
      continue;
    }
    
    if (line.find("Scale factors:") != std::string::npos)
    {
      scaleFactorsFollow = true;
      continue;
    }
    
    if (line.find("#Scale factor sets") != std::string::npos)
    {
      scaleFactorsFollow = false;
      continue;
    }
    
    if (nodesFollow)
    {
      trim(line);
      VLOG(2) << "parse line with nodes: [" << line << "]";
      
      while(!line.empty())
      {
        node_no_t nodeGlobalNo = atoi(line.c_str())-1;
        //VLOG(2) << "       node: " << nodeGlobalNo;
        elements_[elementNo].nodeGlobalNo.push_back(nodeGlobalNo);
        
        // proceed to next number
        if (line.find(" ") != std::string::npos)
        {
          extractUntil(line, " ");
          trim(line);
        }
        else break;
      }
      nodesFollow = false;
    }
    
    if (scaleFactorsFollow)
    {
      trim(line);
      while(!line.empty())
      {
        double scaleFactor = atof(line.c_str());
        //VLOG(2) << "       scaleFactor: " << scaleFactor;
        elements_[elementNo].scaleFactors.push_back(scaleFactor);
        // proceed to next number
        if (line.find(" ") != std::string::npos)
        {
          extractUntil(line, " ");
          trim(line);
        }
        else break;
      }
    }
  }
}

ElementToNodeMapping::Element& ElementToNodeMapping::getElement(element_no_t elementGlobalNo)
{
  assert(elementGlobalNo < (int)elements_.size());
  return elements_[elementGlobalNo];
}

void ElementToNodeMapping::outputElementExelemFile(std::ostream &file, element_no_t elementGlobalNo)
{
  assert(elementGlobalNo < (int)elements_.size());
  
  file << " Element: " << elementGlobalNo+1 << " 0 0" << std::endl 
    << " Nodes:" << std::endl;
    
  // output global node numbers
  for (unsigned int nodeNo = 0; nodeNo < elements_[elementGlobalNo].nodeGlobalNo.size(); nodeNo++)
  {
    file << " " << elements_[elementGlobalNo].nodeGlobalNo[nodeNo]+1;
  }
  file << std::endl 
    << " Scale factors:" << std::endl;
    
  // output scale factors
  StringUtility::outputValuesBlock(file, elements_[elementGlobalNo].scaleFactors, 5);
}
  
};
