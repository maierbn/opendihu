#pragma once

#include <Python.h>  // has to be the first included header

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
    std::vector<int> nodeGlobalNo;      //< the global node numbers of the nodes of the element. E.g. for 3D Hermite elements this vector will contain 8 entries.
    std::vector<double> scaleFactors;   //< the scale factors of the element. there is one scale factor per dof. E.g. for 3D Hermite elements this vector will contain 64 entries (8 dofs per node, 8 nodes).
  };

  //! resize internal representation variable to number of elements
  void setNumberElements(element_no_t nElements);

  //! parse a part of the exelem file that describes a single element
  void parseElementFromExelemFile(std::string content);

  //! return the node numbers and scale factors of the element
  Element &getElement(element_no_t elementNo);

  //! output a single element to exelem file stream
  void outputElementExelem(std::ostream &file, element_no_t elementGlobalNo);

private:

  std::vector<Element> elements_;   //< for global element no the nodes and scale factors
};

} // namespace
