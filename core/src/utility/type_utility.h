#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

namespace TypeUtility
{

/** provides a type trait that determines if the type is a std::vector. If so it defines value as true, if not it defines value as false.
 *  This is analogous to e.g. std::is_array .
 */
template<typename Type>
struct isVector : std::false_type {};
 
template<typename Type, typename Allocator>
struct isVector<std::vector<Type,Allocator>> : std::true_type {};
 
/** Determine if the type is a tuple
 */
template<typename... Types>
struct isTuple : std::false_type {};
 
template<typename... Types>
struct isTuple<std::tuple<Types...>> : std::true_type {};
 
};
