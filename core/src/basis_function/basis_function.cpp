#include "basis_function/basis_function.h"

#include "basis_function/lagrange.h"
#include "basis_function/hermite.h"
#include "utility/string_utility.h"

namespace BasisFunction
{

// linear Lagrange
template <>
std::string getBasisRepresentationString<1,Lagrange<1>>()
{
  return StringUtility::multiply<1>("l.Lagrange");
}

template <>
std::string getBasisRepresentationString<2,Lagrange<1>>()
{
  return StringUtility::multiply<2>("l.Lagrange");
}

template <>
std::string getBasisRepresentationString<3,Lagrange<1>>()
{
  return StringUtility::multiply<3>("l.Lagrange");
}

// quadratic Lagrange
template <>
std::string getBasisRepresentationString<1,Lagrange<2>>()
{
  return StringUtility::multiply<1>("q.Lagrange");
}

template <>
std::string getBasisRepresentationString<2,Lagrange<2>>()
{
  return StringUtility::multiply<2>("q.Lagrange");
}

template <>
std::string getBasisRepresentationString<3,Lagrange<2>>()
{
  return StringUtility::multiply<3>("q.Lagrange");
}

// cubic Hermite
template <>
std::string getBasisRepresentationString<1,Hermite>()
{
  return StringUtility::multiply<1>("c.Hermite");
}

template <>
std::string getBasisRepresentationString<2,Hermite>()
{
  return StringUtility::multiply<2>("c.Hermite");
}

template <>
std::string getBasisRepresentationString<3,Hermite>()
{
  return StringUtility::multiply<3>("c.Hermite");
}

}   // namespace