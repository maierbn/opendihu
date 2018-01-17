#pragma once

#include <petscmat.h>
#include <iostream>
#include <memory>
#include <vector>

namespace FieldVariable
{

class ExfileElementRepresentation
{
public:
  
  struct Node
  {
    std::vector<int> valueIndices;
    std::vector<int> scaleFactorIndices;
  };
  
  //! parse current component's exfile representation from file contents
  void parseFromExelemFile(std::string content);
  
  //! compare operator
  bool operator==(const ExfileElementRepresentation &rhs) const;
  
  //! get the information stored for a node
  Node &getNode(int nodeIndex);
  
  //! output the exelem file header to the stream
  void outputHeaderExelemFile(std::ofstream &file);
  
private:
  
  std::vector<Node> node_;   ///< for the nodes that make up an element their value indices and scale factor indices as states in the exelem file
};

};  // namespace