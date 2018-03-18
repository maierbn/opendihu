#include "field_variable/unstructured/exfile_representation.h"

#include <cassert>

#include "easylogging++.h"
#include "utility/string_utility.h"

namespace FieldVariable
{

using namespace StringUtility;
  
void ExfileRepresentation::parseHeaderFromExelemFile(std::string content)
{
  VLOG(1) << "ExfileRepresentation::parseHeaderFromExelemFile";
  
  std::stringstream str;
  for(unsigned int i=0; i<representation_.size(); i++)
  {
    str << representation_[i];
  }
  VLOG(1) << str.str();
  
  currentElementRepresentation_ = std::make_shared<ExfileElementRepresentation>();
  currentElementRepresentation_->parseFromExelemFile(content);
  
  VLOG(1) << " size of representation: " << representation_.size();
  
  // replace by earlier representation object if there exists one that is equal to the just created one
  for (std::vector<std::shared_ptr<ExfileElementRepresentation>>::iterator iter = representation_.begin(); iter != representation_.end(); iter++)
  {
    
    VLOG(1) << "." << (*iter) << "-" << currentElementRepresentation_;
    if (*iter != nullptr)
    {
      if (**iter == *currentElementRepresentation_)
      {
        currentElementRepresentation_ = *iter;
        break;
      }
    }
  }
}

void ExfileRepresentation::parseElementFromExelemFile(std::string content)
{
  VLOG(1) << "ExfileRepresentation::parseElementFromExelemFile ";
  std::stringstream str;
  for(unsigned int i=0; i<representation_.size(); i++)
  {
    str << representation_[i];
  }
  VLOG(1) << str.str();
  
  // parse element no
  int elementNo = getNumberAfterString(content, "Element:")-1;
  if ((element_no_t)representation_.size() < elementNo)
    representation_.resize(elementNo);
  
  VLOG(1) << " assign current representation for element " << elementNo;
  assert(currentElementRepresentation_);
  representation_[elementNo] = currentElementRepresentation_;
}

void ExfileRepresentation::setNumberElements(element_no_t nElements)
{
  representation_.resize(nElements);
}

bool ExfileRepresentation::operator==(const ExfileRepresentation& rhs)
{
  //VLOG(1) << " compare exfile representation, sizes: " << representation_.size() << "!=" << rhs.representation_.size();
    
  if (representation_.size() != rhs.representation_.size())
  { 
    VLOG(1) << " exfile representation size mismatch: " << representation_.size() << "!=" << rhs.representation_.size();
    return false;
  }
  
  // compare every exfile representation
  for (size_t i=0; i<representation_.size(); i++)
  {
    assert(representation_[i]);
    assert(rhs.representation_[i]);
    
    if (!(*representation_[i] == *rhs.representation_[i]))
    {
      VLOG(1) << " exfile rep different for element " << i;
      return false;
    }
  }
  return true;
}

std::shared_ptr<ExfileElementRepresentation> &ExfileRepresentation::getExfileElementRepresentation(element_no_t elementNo)
{
  assert(elementNo < (element_no_t)representation_.size());
  return representation_[elementNo];
}

void ExfileRepresentation::unifyExfileElementRepresentations()
{
  if (representation_.size() < 2)
    return;
  
  for(element_no_t elementGlobalNo1 = 0; elementGlobalNo1 < representation_.size()-1; elementGlobalNo1++)
  {
    for(element_no_t elementGlobalNo2 = elementGlobalNo1+1; elementGlobalNo2 < representation_.size(); elementGlobalNo2++)
    {
      if (*representation_[elementGlobalNo1] == *representation_[elementGlobalNo2])
      {
        representation_[elementGlobalNo2] = representation_[elementGlobalNo1];
      }
    }
  }
}

bool ExfileRepresentation::haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
{
  assert(element1 < (element_no_t)representation_.size());
  assert(element2 < (element_no_t)representation_.size());
  
  VLOG(2) << "exfileRpers: " << representation_[element1] << "," << representation_[element2];
  return representation_[element1] == representation_[element2];
}

void ExfileRepresentation::output(std::ostream &stream) const
{
  stream << "ExfileRepresentation: size: " << representation_.size() << ": ";
  for(auto &representation : representation_)
  {
    stream << representation << "-";
  }
  stream << std::endl;
  int i=0;
  for(auto &representation : representation_)
  {
    if (representation == nullptr)
      stream << "el " <<i << "  null" << std::endl;
    else
      stream << "el " <<i << "  " << *representation << std::endl;
    i++;
  }
}

std::ostream &operator<<(std::ostream &stream, const ExfileRepresentation &rhs)
{
  rhs.output(stream);
  return stream; 
}

};