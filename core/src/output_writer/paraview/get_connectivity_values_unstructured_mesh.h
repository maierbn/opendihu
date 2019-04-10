#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/function_space.h"

namespace OutputWriter
{

template<typename FunctionSpaceType>
class GetConnectivityValuesUnstructuredMesh
{
public:
  static void get(std::shared_ptr<FunctionSpaceType> functionSpace, std::vector<int> &connectivityValues);
};

template<int D, typename BasisFunctionType>
class GetConnectivityValuesUnstructuredMesh<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>
{
public:

  typedef FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType> FunctionSpaceType;

  static void get(std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>> functionSpace, std::vector<int> &connectivityValues);
};

}  // namespace OutputWriter

#include "output_writer/paraview/get_connectivity_values_unstructured_mesh.tpp"
