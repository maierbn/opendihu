#pragma once

#include <Python.h>  // has to be the first included header
#include <array>

#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "control/types.h"

namespace Mesh
{

/** 
 * A regular mesh with given number of elements in each coordinate direction. The elements
 * have a all the same length/quadratic/cubic size (=mesh width), however the value of the mesh width is not part of this class,
 * but stored by the geometry field.
 * This mesh cannot deform, i.e. it cannot be used for structural mechanics with deformations.
 */
template<int D>
class StructuredRegularFixedOfDimension : public Structured<D>  // StructuredRegularFixed structured_regular_fixed.h
{
public:
  //! constructor from base class
  using Structured<D>::Structured;
 
  //! get mesh width
  virtual double meshWidth() const = 0;  // defined in field_variable/structured/05_field_variable_data_structured_regular_fixed.h
  
private:
 
};

}  // namespace
