#pragma once

#include <Python.h>   // this has to be the first included header
#include <iostream>

namespace BasisFunction
{

class BasisFunction
{
public:
  virtual ~BasisFunction() {}
private:
};

//! return a basis function type string as used in exfiles, e.g. c.Hermite*c.Hermite
template <int D, typename BasisFunctionType>
std::string getBasisRepresentationString();

}  // namespace