#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <iostream>
#include <memory>

#include "field_variable/exfile_element_representation.h"
#include "control/types.h"

namespace FieldVariable
{

class ExfileRepresentation
{
public:
  
  //! resize internal representation variable to number of elements
  void setNumberElements(int nElements);
 
  //! parse current component's exfile representation from file contents
  void parseHeaderFromExelemFile(std::string content);
  
  //! parse a part of the exelem file that describes a single element
  void parseElementFromExelemFile(std::string content);
  
  //! comparison operator
  bool operator==(const ExfileRepresentation &rhs);
  
  //! get the exfile element representation for an element
  std::shared_ptr<ExfileElementRepresentation> getExfileElementRepresentation(int elementNo);

  //! determine if the two elements have the same exfile representation
  bool haveSameExfileRepresentation(element_idx_t element1, element_idx_t element2);
  
  //! output string representation
  void output(std::ostream &stream) const;
  
private:
  std::shared_ptr<ExfileElementRepresentation> currentElementRepresentation_;   ///< the most recent element indexing
 
  std::vector<std::shared_ptr<ExfileElementRepresentation>> representation_;   ///< for every element the exfile representation, i.e. the indices to interpret value blocks in exelem files
};

// output operator
std::ostream &operator<<(std::ostream &stream, const ExfileRepresentation &rhs);

};  // namespace
