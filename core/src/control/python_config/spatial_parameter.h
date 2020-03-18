#pragma once

#include <Python.h>  // has to be the first included header

/** This represents a parameter value given in the python settings file that varies over the computational domain.
 *  An example is a diffusion tensor which is different at every element.
 *
 *  ValueType is the type of the parameter. E.g. double or a matrix for a tensor field.
 *  The parameter is given in the settings file, potentially as a list. If it is a list, the entries correspond to the local or global elements (depending on inputMeshIsGlobal).
 *  Or for a composite mesh, the entries can correspond to the submeshes. Examples:
 *
 *  "inputMeshIsGlobal": True,          # this means all element numbers are the global numbers, False is also possible, then the given values only correspond to the local subdomain
 *  "prefactor": 10,                    # set the prefactor to 10 in all elements
 *  "prefactor": [1,1,1,2,2,2,3,3,3],   # set the prefactor to 1 in elements 0-2, 2 in elements 3-5 and 3 in elements 6-8
 *  "prefactor": {0:5, -1:5},           # set the prefactor to 5 at element 0 and in the last element and to the default value (1) in all other elements
 *  "prefactor": [1,2],                 # when using a submesh, set prefactor to 1 in all elements of the first submesh and to 2 in all elements of the second submesh
 *
 *  Note, if there is a list, the interpretation whether it should be for elements or submeshes depends on the number of entries. The number of entries has to match
 *  either the global number of elements (inputMeshIsGlobal:True), local number of elements (inputMeshIsGlobal:False) or the the number of submeshes.
 *  In case of submeshes, it is still possible to give values for every element of the total composite mesh.
 *
 */
template <typename FunctionSpaceType, typename ValueType>
class SpatialParameter
{
public:

  //! initialize the spatial parameter from the settings, at the given parameter name "keyString"
  void initialize(PythonConfig specificSettings, std::string keyString, ValueType defaultValue, std::shared_ptr<FunctionSpaceType> functionSpace);

  //! get the value of the parameter in the specified element
  void getValue(element_no_t elementNoLocal, ValueType &value) const;

  //! the value of the parameter in the given element
  //ValueType value(element_no_t elementNoLocal) const;

  //! the value of the parameter in the given element
  const ValueType &value(element_no_t elementNoLocal) const;

protected:

  std::vector<int> valueIndices_;        //< for each elementNoLocal the index into values of the parameter
  std::vector<ValueType> values_;   //< all stored parameter values, to which element they belong is given by valueNo

  bool inputMeshIsGlobal_;          //< if the input is given for all global element numbers
};


/** helper class that sets the indices to the values */
template<typename FunctionSpaceType, typename ValueType>
class SetValueNos
{
public:

  //! set indices in valueIndices_ according to the number of values that were present in the config
  static void set(std::shared_ptr<FunctionSpaceType> functionSpace, bool inputMeshIsGlobal, std::string path,
                  std::vector<ValueType> &values, std::vector<int> &valueIndices);
};

/** implementation for composite meshes */
template<int D, typename BasisFunctionType, typename ValueType>
class SetValueNos<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,ValueType>
{
public:

  //! set indices in valueIndices_ according to the number of values that were present in the config
  static void set(std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>> functionSpace, bool inputMeshIsGlobal, std::string path,
                  std::vector<ValueType> &values, std::vector<int> &valueIndices);
};

#include "control/python_config/spatial_parameter.tpp"
