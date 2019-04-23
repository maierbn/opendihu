#pragma once

#include <Python.h>  // has to be the first included header

namespace FunctionSpace
{

// define generic function space without logical real world mesh presententation, that can be used for generic field variables.
// For example for MOR the reduced vectors do not live on any mesh, but they need a function space to be defined and such that output writers work.
typedef FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> Generic;

}  // namespace
