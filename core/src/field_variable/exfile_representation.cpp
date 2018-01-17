#include "field_variable/exfile_representation.h"

#include <cassert>

#include "utility/string_utility.h"

namespace FieldVariable
{

using namespace StringUtility;
  
void ExfileRepresentation::parseHeaderFromExelemFile(std::string content)
{
  currentElementRepresentation_ = std::make_shared<ExfileElementRepresentation>();
  currentElementRepresentation_->parseFromExelemFile(content);
  
  // replace by earlier representation object if there exists one that is equal to the just created one
  for (auto &representation : representation_)
  {
    if (*representation == *currentElementRepresentation_)
    {
      currentElementRepresentation_ = representation;
      break;
    }
  }
}

void ExfileRepresentation::parseElementFromExelemFile(std::string content)
{
  int elementNo = getNumberAfterString(content, "Element:");
  if ((int)representation_.size() < elementNo)
    representation_.resize(elementNo);
  representation_[elementNo-1] = currentElementRepresentation_;
}

void ExfileRepresentation::setNumberElements(int nElements)
{
  representation_.resize(nElements);
}

bool ExfileRepresentation::operator==(const ExfileRepresentation& rhs)
{
  if (representation_.size() != rhs.representation_.size())
    return false;
  for (size_t i=0; i<representation_.size(); i++)
  {
    if (!(*representation_[i] == *rhs.representation_[i]))
      return false;
  }
  return true;
}

std::shared_ptr<ExfileElementRepresentation> ExfileRepresentation::getExfileElementRepresentation(int elementNo)
{
  if (elementNo >= (int)representation_.size())
    return nullptr;
  return representation_[elementNo];
}

bool ExfileRepresentation::haveSameExfileRepresentation(element_idx_t element1, element_idx_t element2)
{
  assert(element1 < representation_.size());
  assert(element2 < representation_.size());
  return representation_[element1] == representation_[element2];
}

};