#include "function_space/mixed_function_space.h"

#include <Python.h>  // has to be the first included header

namespace FunctionSpace
{

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType>
Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>::
Mixed(PyObject *specificSettings) : Mesh::Mesh(specificSettings),
  lowOrderFunctionSpace_(std::make_shared<LowOrderFunctionSpaceType>(specificSettings, true)),
  highOrderFunctionSpace_(std::make_shared<HighOrderFunctionSpaceType>(specificSettings, false))
{
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType>
void Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>::
initialize()
{
  lowOrderFunctionSpace_->initialize();
  highOrderFunctionSpace_->initialize();
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType>
int Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>::
dimension() const
{
  return highOrderFunctionSpace_->dimension();
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType>
constexpr int Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>::
dim()
{
  return HighOrderFunctionSpaceType::Mesh::dim();
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType>
node_no_t Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>::
nNodesLocalWithGhosts() const
{
  return highOrderFunctionSpace_->nNodesLocalWithGhosts();
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType>
std::shared_ptr<LowOrderFunctionSpaceType> Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>::
lowOrderFunctionSpace()
{
  return lowOrderFunctionSpace_;
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType>
std::shared_ptr<HighOrderFunctionSpaceType> Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>::
highOrderFunctionSpace()
{
  return highOrderFunctionSpace_;
}

};  //namespace