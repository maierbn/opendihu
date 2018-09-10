#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <memory>

namespace BoundaryConditions
{

//! parse config and extract boundary conditions specified under key "dirichletBoundaryConditions", T can be std::array<double,nComponents>
template<typename FunctionSpaceType, typename T>
void parseBoundaryConditions(PyObject *settings, std::shared_ptr<FunctionSpaceType> functionSpace,
                             std::vector<std::pair<int,T>> &boundaryConditions);

};  // namespace

#include "utility/boundary_conditions.tpp"
