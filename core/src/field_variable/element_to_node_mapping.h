#pragma once

#include <petscmat.h>
#include <iostream>
#include <memory>
#include <vector>
#include <fstream>

#include "control/types.h"

namespace FieldVariable
{

class ElementToNodeMapping
{
public:
 
  /** global node numbers and scale factors for each element
   */
  struct Element 
  {
    std::vector<int> globalNodeNo;
    std::vector<double> scaleFactors;
  };
 
  //! resize internal representation variable to number of elements
  void setNumberElements(int nElements);
 
  //! parse a part of the exelem file that describes a single element
  void parseElementFromExelemFile(std::string content);

  //! return the node numbers and scale factors of the element
  Element &getElement(int elementNo);
  
  //! output a single element to exelem file stream
  void outputElementExelemFile(std::ofstream &file, element_idx_t elementGlobalNo);
  
private:
 
  std::vector<Element> elements_;
};

};  // namespace