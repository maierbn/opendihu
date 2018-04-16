#pragma once

#include <Python.h>  // has to be the first included header

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
    std::vector<int> valueIndices;        ///< the indices of the dof values of this node in the exnode file node values (sub-)block (for the particular field variable/component) of the node. If there are not multiple versions, this is simply 0,1,...,ndofs-1. If there are e.g. 2 versions and 8 dofs per node, this can be 0,1,...,7 if the elements uses the 1st version, or 8,...,15 if the element uses the second version. Note, that the real index of the dofs inside the values block may be different when this is not the first component of the block.
    std::vector<int> scaleFactorIndices;   ///< the indices of all scale factor entries for this node in the exelem element scale factors block. Thus this is kind of a node to element block mapping. 
  };
  
  //! parse current component's exfile representation from file contents
  void parseFromExelemFile(std::string content);
 
  //! resize the internal node_ vector to number of nodes
  void setNumberNodes(int nNodes);
  
  //! comparison operator
  bool operator==(const ExfileElementRepresentation &rhs) const;
  
  //! get the information stored for a node, identified by local index
  Node &getNode(int nodeIndex);
  
  //! output the exelem file header to the stream
  void outputHeaderExelem(std::ostream &file);
  
  //! output string representation
  void output(std::ostream &stream) const;
  
private:
  
  std::vector<Node> node_;   ///< for the nodes that make up an element their value indices and scale factor indices as stated in the exelem file
};

// output operator
std::ostream &operator<<(std::ostream &stream, const ExfileElementRepresentation &rhs);

};  // namespace
