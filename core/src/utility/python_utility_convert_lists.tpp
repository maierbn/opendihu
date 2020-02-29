#include "utility/python_utility.h"

#include <Python.h>
#include "easylogging++.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>

#include "utility/vector_operators.h"
#include "control/types.h"
#include "utility/matrix.h"


//partial specialization for std::array
template<typename ValueType, int nComponents>
struct PythonUtility::convertFromPython<std::array<ValueType,nComponents>>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static std::array<ValueType,nComponents> get(PyObject *object, std::array<ValueType,nComponents> defaultValue)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    std::array<ValueType,nComponents> result;
    assert(object != nullptr);
    if (PyList_Check(object))
    {
      int i = 0;
      int iEnd = std::min((int)PyList_Size(object), nComponents);

      for (;i < iEnd; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyList_GetItem(object, (Py_ssize_t)i), defaultValue[i]);
      }

      // fill rest of values with default values
      for (;i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }

      if (iEnd < nComponents)
      {
        LOG(WARNING) << "Python list only contains " << iEnd << " values, but " << nComponents << " are required. Filling rest with default values."
          << " Parsed values: " << result;
      }

      return result;
    }
    else if (PyTuple_Check(object))
    {
      int i = 0;
      int iEnd = std::min((int)PyTuple_Size(object), nComponents);

      for (;i < iEnd; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyTuple_GetItem(object, (Py_ssize_t)i), defaultValue[i]);
      }

      // fill rest of values with default values
      for (;i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }

      if (iEnd < nComponents)
      {
        LOG(WARNING) << "Python tuple only contains " << iEnd << " values, but " << nComponents << " are required. Filling rest with default values."
          << " Parsed values: " << result;
      }
      return result;
    }
    else if (PyDict_Check(object))
    {
      PyObject *itemList = PyDict_Items(object);

      for (int itemListIndex = 0; itemListIndex < PyList_Size(itemList); itemListIndex++)
      {
        PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
        PyObject *pyKey = PyTuple_GetItem(tuple, (Py_ssize_t)0);
        PyObject *pyValue = PyTuple_GetItem(tuple, (Py_ssize_t)1);

        int key = convertFromPython<int>::get(pyKey);
        if (key >= 0 && key < nComponents)
        {
          result[key] = convertFromPython<ValueType>::get(pyValue);
        }
      }
      return result;
    }
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object, defaultValue[0]);

      result[0] = valueDouble;
      std::copy(defaultValue.begin()+1, defaultValue.end(), result.begin()+1);

      if (nComponents > 1)
      {
        LOG(WARNING) << "Expected python list with " << nComponents << " entries, but only a single value was found. Filling rest with default values."
          << " Parsed values: " << result;
      }
      return result;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static std::array<ValueType,nComponents> get(PyObject *object)
  {
    std::array<ValueType,nComponents> defaultValue;
    defaultValue.fill(ValueType());
    return convertFromPython<std::array<ValueType,nComponents>>::get(object, defaultValue);
  }
};

//partial specialization for std::array
template<typename ValueType, global_no_t nComponents>
struct PythonUtility::convertFromPython<std::array<ValueType,nComponents>>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static std::array<ValueType,nComponents> get(PyObject *object, std::array<ValueType,nComponents> defaultValue)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    std::array<ValueType,nComponents> result;
    assert(object != nullptr);
    if (PyList_Check(object))
    {
      global_no_t i = 0;
      global_no_t iEnd = std::min((global_no_t)PyList_Size(object), nComponents);

      for (;i < iEnd; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyList_GetItem(object, (Py_ssize_t)i), defaultValue[i]);
      }

      // fill rest of values with default values
      for (;i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }
      if (iEnd < nComponents)
      {
        LOG(WARNING) << "Python list only contains " << iEnd << " values, but " << nComponents << " are required. Filling rest with default values."
          << " Parsed values: " << result;
      }
      return result;
    }
    else if (PyTuple_Check(object))
    {
      int i = 0;
      int iEnd = std::min((global_no_t)PyTuple_Size(object), nComponents);

      for (;i < iEnd; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyTuple_GetItem(object, (Py_ssize_t)i), defaultValue[i]);
      }

      // fill rest of values with default values
      for (;i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }
      if (iEnd < nComponents)
      {
        LOG(WARNING) << "Python tuple only contains " << iEnd << " values, but " << nComponents << " are required. Filling rest with default values."
          << " Parsed values: " << result;
      }
      return result;
    }
    else if (PyDict_Check(object))
    {
      PyObject *itemList = PyDict_Items(object);

      for (int itemListIndex = 0; itemListIndex < PyList_Size(itemList); itemListIndex++)
      {
        PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
        PyObject *pyKey = PyTuple_GetItem(tuple, (Py_ssize_t)0);
        PyObject *pyValue = PyTuple_GetItem(tuple, (Py_ssize_t)1);

        int key = convertFromPython<int>::get(pyKey);
        if (key >= 0 && key < nComponents)
        {
          result[key] = convertFromPython<ValueType>::get(pyValue);
        }
      }
      return result;
    }
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object, defaultValue[0]);

      result[0] = valueDouble;
      std::copy(defaultValue.begin()+1, defaultValue.end(), result.begin()+1);

      if (nComponents > 1)
      {
        LOG(WARNING) << "Expected python list with " << nComponents << " entries, but only a single value was found. Filling rest with default values."
          << " Parsed values: " << result;
      }
      return result;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static std::array<ValueType,nComponents> get(PyObject *object)
  {
    std::array<ValueType,nComponents> defaultValue;
    defaultValue.fill(ValueType());
    return convertFromPython<std::array<ValueType,nComponents>>::get(object, defaultValue);
  }
};

//partial specialization for std::array
template<typename ValueType, unsigned long nComponents>
struct PythonUtility::convertFromPython<std::array<ValueType,nComponents>>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static std::array<ValueType,nComponents> get(PyObject *object, std::array<ValueType,nComponents> defaultValue)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    std::array<ValueType,nComponents> result;
    assert(object != nullptr);
    if (PyList_Check(object))
    {
      unsigned long i = 0;
      unsigned long iEnd = std::min((unsigned long)PyList_Size(object), nComponents);

      for (;i < iEnd; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyList_GetItem(object, (Py_ssize_t)i), (ValueType)defaultValue[i]);
      }

      // fill rest of values with default values
      for (;i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }
      if (iEnd < nComponents)
      {
        LOG(WARNING) << "Python list only contains " << iEnd << " values, but " << nComponents << " are required. Filling rest with default values."
          << " Parsed values: " << result;
      }
      return result;
    }
    else if (PyTuple_Check(object))
    {
      int i = 0;
      int iEnd = std::min((unsigned long)PyTuple_Size(object), nComponents);

      for (;i < iEnd; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyTuple_GetItem(object, (Py_ssize_t)i), defaultValue[i]);
      }

      // fill rest of values with default values
      for (;i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }
      if (iEnd < nComponents)
      {
        LOG(WARNING) << "Python tuple only contains " << iEnd << " values, but " << nComponents << " are required. Filling rest with default values."
          << " Parsed values: " << result;
      }
      return result;
    }
    else if (PyDict_Check(object))
    {
      PyObject *itemList = PyDict_Items(object);

      for (int itemListIndex = 0; itemListIndex < PyList_Size(itemList); itemListIndex++)
      {
        PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
        PyObject *pyKey = PyTuple_GetItem(tuple, (Py_ssize_t)0);
        PyObject *pyValue = PyTuple_GetItem(tuple, (Py_ssize_t)1);

        int key = convertFromPython<int>::get(pyKey);
        if (key >= 0 && key < nComponents)
        {
          result[key] = convertFromPython<ValueType>::get(pyValue);
        }
      }
      return result;
    }
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object, defaultValue[0]);

      result[0] = valueDouble;
      std::copy(defaultValue.begin()+1, defaultValue.end(), result.begin()+1);

      if (nComponents > 1)
      {
        LOG(WARNING) << "Expected python list with " << nComponents << " entries, but only a single value was found. Filling rest with default values."
          << " Parsed values: " << result;
      }
      return result;
    }
#ifndef __PGI
    return defaultValue;
#endif
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static std::array<ValueType,nComponents> get(PyObject *object)
  {
    std::array<ValueType,nComponents> defaultValue;
    defaultValue.fill(ValueType());
    return convertFromPython<std::array<ValueType,nComponents>>::get(object, defaultValue);
  }
};

//partial specialization for MathUtility::Matrix
template<int nRows, int nColumns>
struct PythonUtility::convertFromPython<MathUtility::Matrix<nRows,nColumns>>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static MathUtility::Matrix<nRows,nColumns> get(PyObject *object, MathUtility::Matrix<nRows,nColumns> defaultValue)
  {
    return PythonUtility::convertFromPython<std::array<double, nRows*nColumns>>::get(object, defaultValue);
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static MathUtility::Matrix<nRows,nColumns> get(PyObject *object)
  {
    return convertFromPython<std::array<double, nRows*nColumns>>::get(object);
  }
};

//partial specialization for std::vector
template<typename ValueType>
struct PythonUtility::convertFromPython<std::vector<ValueType>>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static std::vector<ValueType> get(PyObject *object, std::vector<ValueType> defaultValue)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    std::vector<ValueType> result;
    assert(object != nullptr);
    if (PyList_Check(object))
    {
      int nEntries = (int)PyList_Size(object);
      result.resize(nEntries);
      LOG(DEBUG) << "nEntries: " << nEntries;

      for (int i = 0; i < nEntries; i++)
      {
        ValueType defaultValueItem;
        if (defaultValue.size() >= i)
        {
          defaultValueItem = defaultValue[i];
        }
        else if (defaultValue.size() >= 1)
        {
          defaultValueItem = defaultValue[0];
        }
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyList_GetItem(object, (Py_ssize_t)i), defaultValueItem);
      }
      return result;
    }
    else if (PyTuple_Check(object))
    {
      int nEntries = PyTuple_Size(object);
      result.resize(nEntries + 1);

      for (int i = 0; i < nEntries; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyTuple_GetItem(object, (Py_ssize_t)i), defaultValue[i]);
      }
      return result;
    }
    else if (PyDict_Check(object))
    {
      PyObject *itemList = PyDict_Items(object);

      for (int itemListIndex = 0; itemListIndex < PyList_Size(itemList); itemListIndex++)
      {
        PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
        PyObject *pyKey = PyTuple_GetItem(tuple, (Py_ssize_t)0);
        PyObject *pyValue = PyTuple_GetItem(tuple, (Py_ssize_t)1);

        int key = convertFromPython<int>::get(pyKey);
        if (key >= 0)
        {
          if (key >= result.size())
            result.resize(key + 1);
          result[key] = convertFromPython<ValueType>::get(pyValue);
        }
      }
      return result;
    }
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object, defaultValue[0]);

      result.resize(1);
      result[0] = valueDouble;

      return result;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking
  static std::vector<ValueType> get(PyObject *object, ValueType defaultValue)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    std::vector<ValueType> result;
    assert(object != nullptr);
    if (PyList_Check(object))
    {
      int nEntries = (int)PyList_Size(object);
      result.resize(nEntries);

      for (int i = 0; i < nEntries; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyList_GetItem(object, (Py_ssize_t)i));
      }
      return result;
    }
    else if (PyTuple_Check(object))
    {
      int nEntries = PyTuple_Size(object);
      result.resize(nEntries + 1);

      for (int i = 0; i < nEntries; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyTuple_GetItem(object, (Py_ssize_t)i));
      }
      return result;
    }
    else if (PyDict_Check(object))
    {
      PyObject *itemList = PyDict_Items(object);

      for (int itemListIndex = 0; itemListIndex < PyList_Size(itemList); itemListIndex++)
      {
        PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
        PyObject *pyKey = PyTuple_GetItem(tuple, (Py_ssize_t)0);
        PyObject *pyValue = PyTuple_GetItem(tuple, (Py_ssize_t)1);

        int key = convertFromPython<int>::get(pyKey);

        if (key >= 0)
        {
          if (key >= result.size())
            result.resize(key + 1, defaultValue);
          result[key] = convertFromPython<ValueType>::get(pyValue);
        }
      }
      return result;
    }
    else if (object == Py_None)
    {
      // object is None, this means empty list
      return result;
    }
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object);

      result.resize(1);
      result[0] = valueDouble;

      return result;
    }

#ifndef __PGI
    return std::vector<ValueType>();
#endif
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static std::vector<ValueType> get(PyObject *object)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    std::vector<ValueType> result;
    assert(object != nullptr);
    if (PyList_Check(object))
    {
      int nEntries = (int)PyList_Size(object);
      result.resize(nEntries);

      for (int i = 0; i < nEntries; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyList_GetItem(object, (Py_ssize_t)i));
      }
      return result;
    }
    else if (PyTuple_Check(object))
    {
      int nEntries = PyTuple_Size(object);
      result.resize(nEntries + 1);

      for (int i = 0; i < nEntries; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyTuple_GetItem(object, (Py_ssize_t)i));
      }
      return result;
    }
    else if (PyDict_Check(object))
    {
      PyObject *itemList = PyDict_Items(object);

      for (int itemListIndex = 0; itemListIndex < PyList_Size(itemList); itemListIndex++)
      {
        PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
        PyObject *pyKey = PyTuple_GetItem(tuple, (Py_ssize_t)0);
        PyObject *pyValue = PyTuple_GetItem(tuple, (Py_ssize_t)1);

        int key = convertFromPython<int>::get(pyKey);

        if (key >= 0)
        {
          if (key >= result.size())
            result.resize(key + 1, convertFromPython<ValueType>::get(Py_None));
          result[key] = convertFromPython<ValueType>::get(pyValue);
        }
      }
      return result;
    }
    else if (object == Py_None)
    {
      // object is None, this means empty list
      return result;
    }
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object);

      result.resize(1);
      result[0] = valueDouble;

      return result;
    }

#ifndef __PGI
    return std::vector<ValueType>();
#endif
  }
};

//partial specialization for std::vector<std::pair<>>, this can convert tuples to [<key,value>,...] vectors
template<typename KeyType,typename ValueType>
struct PythonUtility::convertFromPython<std::vector<std::pair<KeyType,ValueType>>>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static std::vector<std::pair<KeyType,ValueType>> get(PyObject *object)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    std::vector<std::pair<KeyType,ValueType>> result;
    assert(object != nullptr);

    if (PyDict_Check(object))
    {
      PyObject *itemList = PyDict_Items(object);

      for (int itemListIndex = 0; itemListIndex < PyList_Size(itemList); itemListIndex++)
      {
        PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
        PyObject *pyKey = PyTuple_GetItem(tuple, (Py_ssize_t)0);
        PyObject *pyValue = PyTuple_GetItem(tuple, (Py_ssize_t)1);

        result.push_back(std::make_pair<KeyType,ValueType>(
          convertFromPython<KeyType>::get(pyKey),
          convertFromPython<ValueType>::get(pyValue)
        ));
      }
      return result;
    }
    else if (PyList_Check(object))
    {
      int nEntries = (int)PyList_Size(object);

      std::vector<std::pair<KeyType,ValueType>> result(nEntries);

      for (int i = 0; i < nEntries; i++)
      {
        result[i] = PythonUtility::convertFromPython<std::pair<KeyType,ValueType>>::get(PyList_GetItem(object, (Py_ssize_t)i));
      }
      return result;
    }
    else if (PyTuple_Check(object))
    {
      int nEntries = PyTuple_Size(object);

      std::vector<std::pair<KeyType,ValueType>> result(nEntries);

      for (int i = 0; i < nEntries; i++)
      {
        result[i] = PythonUtility::convertFromPython<std::pair<KeyType,ValueType>>::get(PyTuple_GetItem(object, (Py_ssize_t)i));
      }
      return result;
    }
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object);

      std::vector<std::pair<KeyType,ValueType>> result(1);
      result[0] = std::pair<KeyType,ValueType>(KeyType(),valueDouble);

      return result;
    }

#ifndef __PGI
    return std::vector<std::pair<KeyType,ValueType>>();
#endif
  }
};

//partial specialization for std::pair
template<typename ValueType1,typename ValueType2>
struct PythonUtility::convertFromPython<std::pair<ValueType1,ValueType2>>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static std::pair<ValueType1,ValueType2> get(PyObject *object, std::pair<ValueType1,ValueType2> defaultValue)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    std::pair<ValueType1,ValueType2> result;
    assert(object != nullptr);
    if (PyTuple_Check(object))
    {
      int nEntries = (int)PyTuple_Size(object);
      if (nEntries != 2)
      {
        LOG(WARNING) << "Converting python tuple to pair, but the tuple has " << nEntries << " != 2 entries.";
      }

      result.first = PythonUtility::convertFromPython<ValueType1>::get(PyTuple_GetItem(object, (Py_ssize_t)0), defaultValue.first);
      result.second = PythonUtility::convertFromPython<ValueType2>::get(PyTuple_GetItem(object, (Py_ssize_t)1), defaultValue.second);

      return result;
    }
    else
    {
      LOG(ERROR) << "Could not convert python tuple to pair.";
    }

    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static std::pair<ValueType1,ValueType2> get(PyObject *object)
  {
    std::pair<ValueType1,ValueType2> defaultValue;
    return PythonUtility::convertFromPython<std::pair<ValueType1,ValueType2>>::get(object, defaultValue);
  }
};

//partial specialization for std::tuple
template <size_t index, typename... ValueTypes>
typename std::enable_if<index == sizeof...(ValueTypes)>::type
ExtractTuple(PyObject *object, int nEntries, std::tuple<ValueTypes...> &result)
{}

template <size_t index, typename... ValueTypes>
typename std::enable_if<(index < sizeof...(ValueTypes))>::type
ExtractTuple(PyObject *object, int nEntries, std::tuple<ValueTypes...> &result)
{
  typedef typename std::tuple_element<index,std::tuple<ValueTypes...>>::type Type;
  if (index < nEntries)
  {
    // extract value from python config
    std::get<index>(result) = PythonUtility::convertFromPython<Type>::get(
      PyTuple_GetItem(object, (Py_ssize_t)index));
  }
  else
  {
    // if the current index is higher than the number of entries in the python config, default values are used
    std::get<index>(result) = Type();
  }

  ExtractTuple<index+1,ValueTypes...>(object, nEntries, result);
}

template<typename... ValueTypes>
struct PythonUtility::convertFromPython<std::tuple<ValueTypes...>>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static std::tuple<ValueTypes...> get(PyObject *object)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    std::tuple<ValueTypes...> result;
    assert(object != nullptr);
    if (PyTuple_Check(object))
    {
      int nEntries = (int)PyTuple_Size(object);

      int nEntriesTuple = std::tuple_size<std::tuple<ValueTypes...>>::value;
      if (nEntriesTuple < nEntries)
      {
        LOG(WARNING) << "Python tuple has " << nEntries << " entries, expected are only " << nEntriesTuple << ".";
      }
      else if (nEntriesTuple > nEntries)
      {
        LOG(WARNING) << "Python tuple has only " << nEntries << " entries, but " << nEntriesTuple << " are required, using 0 for the other values.";
      }

      ExtractTuple<0, ValueTypes...>(object, nEntries, result);
    }
    else
    {
      LOG(ERROR) << "Could not convert python tuple to pair.";
    }

    return result;
  }
};

// ------------- convertToPython ----------------

template<typename ValueType>
struct PythonUtility::convertToPython<std::vector<ValueType>>
{
  //! convert a type to a python object
  static PyObject *get(std::vector<ValueType> &value)
  {
    PyObject *result = PyList_New(value.size());
    for (int i = 0; i < value.size(); i++)
    {
      PyObject *item = PythonUtility::convertToPython<ValueType>::get(value[i]);
      PyList_SET_ITEM(result, (Py_ssize_t)i, item);
    }
    return result;
  }
};

template<typename ValueType, int nComponents>
struct PythonUtility::convertToPython<std::array<ValueType,nComponents>>
{
  //! convert a type to a python object
  static PyObject *get(std::array<ValueType,nComponents> &value)
  {
    PyObject *result = PyTuple_New(nComponents);
    for (int i = 0; i < nComponents; i++)
    {
      PyObject *item = PythonUtility::convertToPython<ValueType>::get(value[i]);
      PyTuple_SET_ITEM(result, (Py_ssize_t)i, item);
    }
    return result;
  }
};

template<typename ValueType, global_no_t nComponents>
struct PythonUtility::convertToPython<std::array<ValueType,nComponents>>
{
  //! convert a type to a python object
  static PyObject *get(std::array<ValueType,nComponents> &value)
  {
    PyObject *result = PyTuple_New(nComponents);
    for (int i = 0; i < nComponents; i++)
    {
      PyObject *item = PythonUtility::convertToPython<ValueType>::get(value[i]);
      PyTuple_SET_ITEM(result, (Py_ssize_t)i, item);
    }
    return result;
  }
};

template<typename ValueType, unsigned long nComponents>
struct PythonUtility::convertToPython<std::array<ValueType,nComponents>>
{
  //! convert a type to a python object
  static PyObject *get(std::array<ValueType,nComponents> &value)
  {
    PyObject *result = PyTuple_New(nComponents);
    for (int i = 0; i < nComponents; i++)
    {
      PyObject *item = PythonUtility::convertToPython<ValueType>::get(value[i]);
      PyTuple_SET_ITEM(result, (Py_ssize_t)i, item);
    }
    return result;
  }
};
