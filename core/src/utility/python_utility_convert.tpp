#include "utility/python_utility.h"

#include <Python.h>
#include "easylogging++.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>

#include "control/use_numpy.h"
#include "utility/vector_operators.h"
#include "control/types.h"

//partial specialization for std::array
template<typename ValueType, int nComponents>
struct PythonUtility::convertFromPython<std::array<ValueType,nComponents>>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static std::array<ValueType,nComponents> get(PyObject *object, std::array<ValueType,nComponents> defaultValue)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    std::array<ValueType,nComponents> result;
    if (PyList_Check(object))
    {
      int i = 0;
      int iEnd = std::min((int)PyList_Size(object), nComponents);

      for(;i < iEnd; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyList_GetItem(object, (Py_ssize_t)i), defaultValue[i]);
      }

      // fill rest of values with default values
      for(;i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }
      return result;
    }
    /*
  #ifdef HAVE_NUMPYC
    else if (PyArray_Check(object))
    {
      LOG(DEBUG) << "object is a pyarray ";

      const PyArrayObject *arrayObject = (PyArrayObject*)object;

      int nElementsInNumpyArray = PyArray_Size(object);
      int nElementsToCopy = std::min(nComponents, nElementsInNumpyArray);

      int typenumber = PyArray_TYPE(arrayObject);  // get the type of the contained data in the numpy array, e.g. NPY_DOUBLE

      if (sizeof(ValueType) == PyArray_ITEMSIZE(arrayObject))
      {
        PyArray_DescrFromType(typenumber)->f->copyswapn(
        result.data(), 1, PyArray_DATA((PyArrayObject*)arrayObject), 1, nElementsToCopy, 0, NULL);
      }
      else
      {
        LOG(ERROR) << "Could not convert numpy array with itemsize " << PyArray_ITEMSIZE(arrayObject) << " bytes to type " << typeid(ValueType).name() << " with size " << sizeof(ValueType) << ".";
        nElementsToCopy = 0;
      }

      // fill rest of values with default values
      for(int i = nElementsToCopy; i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }
      LOG(DEBUG) << "converted to " << nComponents << " entries: " << result;
    }
  #endif
  */
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object, defaultValue[0]);

      result[0] = valueDouble;
      std::copy(defaultValue.begin()+1, defaultValue.end(), result.begin()+1);

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
    PythonUtility::GlobalInterpreterLock lock;

    std::array<ValueType,nComponents> result;
    if (PyList_Check(object))
    {
      global_no_t i = 0;
      global_no_t iEnd = std::min((global_no_t)PyList_Size(object), nComponents);

      for(;i < iEnd; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyList_GetItem(object, (Py_ssize_t)i), defaultValue[i]);
      }

      // fill rest of values with default values
      for(;i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }
      return result;
    }
    /*
  #ifdef HAVE_NUMPYC
    else if (PyArray_Check(object))
    {
      LOG(DEBUG) << "object is a pyarray ";

      const PyArrayObject *arrayObject = (PyArrayObject*)object;

      int nElementsInNumpyArray = PyArray_Size(object);
      int nElementsToCopy = std::min(nComponents, nElementsInNumpyArray);

      int typenumber = PyArray_TYPE(arrayObject);  // get the type of the contained data in the numpy array, e.g. NPY_DOUBLE

      if (sizeof(ValueType) == PyArray_ITEMSIZE(arrayObject))
      {
        PyArray_DescrFromType(typenumber)->f->copyswapn(
        result.data(), 1, PyArray_DATA((PyArrayObject*)arrayObject), 1, nElementsToCopy, 0, NULL);
      }
      else
      {
        LOG(ERROR) << "Could not convert numpy array with itemsize " << PyArray_ITEMSIZE(arrayObject) << " bytes to type " << typeid(ValueType).name() << " with size " << sizeof(ValueType) << ".";
        nElementsToCopy = 0;
      }

      // fill rest of values with default values
      for(int i = nElementsToCopy; i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }
      LOG(DEBUG) << "converted to " << nComponents << " entries: " << result;
    }
  #endif
  */
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object, defaultValue[0]);

      result[0] = valueDouble;
      std::copy(defaultValue.begin()+1, defaultValue.end(), result.begin()+1);

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
    PythonUtility::GlobalInterpreterLock lock;

    std::array<ValueType,nComponents> result;
    if (PyList_Check(object))
    {
      unsigned long i = 0;
      unsigned long iEnd = std::min((unsigned long)PyList_Size(object), nComponents);

      for(;i < iEnd; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyList_GetItem(object, (Py_ssize_t)i), (ValueType)defaultValue[i]);
      }

      // fill rest of values with default values
      for(;i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }
      return result;
    }
    /*
  #ifdef HAVE_NUMPYC
    else if (PyArray_Check(object))
    {
      LOG(DEBUG) << "object is a pyarray ";

      const PyArrayObject *arrayObject = (PyArrayObject*)object;

      int nElementsInNumpyArray = PyArray_Size(object);
      int nElementsToCopy = std::min(nComponents, nElementsInNumpyArray);

      int typenumber = PyArray_TYPE(arrayObject);  // get the type of the contained data in the numpy array, e.g. NPY_DOUBLE

      if (sizeof(ValueType) == PyArray_ITEMSIZE(arrayObject))
      {
        PyArray_DescrFromType(typenumber)->f->copyswapn(
        result.data(), 1, PyArray_DATA((PyArrayObject*)arrayObject), 1, nElementsToCopy, 0, NULL);
      }
      else
      {
        LOG(ERROR) << "Could not convert numpy array with itemsize " << PyArray_ITEMSIZE(arrayObject) << " bytes to type " << typeid(ValueType).name() << " with size " << sizeof(ValueType) << ".";
        nElementsToCopy = 0;
      }

      // fill rest of values with default values
      for(int i = nElementsToCopy; i < nComponents; i++)
      {
        result[i] = defaultValue[i];
      }
      LOG(DEBUG) << "converted to " << nComponents << " entries: " << result;
    }
  #endif
  */
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object, defaultValue[0]);

      result[0] = valueDouble;
      std::copy(defaultValue.begin()+1, defaultValue.end(), result.begin()+1);

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

//partial specialization for std::vector
template<typename ValueType>
struct PythonUtility::convertFromPython<std::vector<ValueType>>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static std::vector<ValueType> get(PyObject *object, std::vector<ValueType> defaultValue)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    std::vector<ValueType> result;
    if (PyList_Check(object))
    {
      int nEntries = (int)PyList_Size(object);
      result.resize(nEntries);

      for(int i = 0; i < nEntries; i++)
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
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object, defaultValue[0]);

      result.resize(1);
      result[0] = valueDouble;

      return result;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static std::vector<ValueType> get(PyObject *object)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    std::vector<ValueType> result;
    if (PyList_Check(object))
    {
      int nEntries = (int)PyList_Size(object);
      result.resize(nEntries);

      for(int i = 0; i < nEntries; i++)
      {
        result[i] = PythonUtility::convertFromPython<ValueType>::get(PyList_GetItem(object, (Py_ssize_t)i));
      }
      return result;
    }
    else
    {
      ValueType valueDouble = PythonUtility::convertFromPython<ValueType>::get(object);

      result.resize(1);
      result[0] = valueDouble;

      return result;
    }
    return std::vector<ValueType>();
  }
};

//partial specialization for int
template<>
struct PythonUtility::convertFromPython<int>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static int get(PyObject *object, int defaultValue)
  {
    if(object == NULL)
      return defaultValue;

    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    if (PyLong_Check(object))
    {
      long valueLong = PyLong_AsLong(object);
      return int(valueLong);
    }
    else if (PyFloat_Check(object))
    {
      double valueDouble = PyFloat_AsDouble(object);

      if (double(int(valueDouble)) != valueDouble)      // if value is not e.g. 2.0
      {
        LOG(WARNING) << "convertFromPython: object is no int: " << object;
      }

      return int(valueDouble);
    }
    else if (PyUnicode_Check(object))
    {
      std::string valueString = pyUnicodeToString(object);
      return atoi(valueString.c_str());
    }
    else
    {
      LOG(WARNING) << "convertFromPython: object is no int: " << object;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static int get(PyObject *object)
  {
    return convertFromPython<int>::get(object, 0);
  }
};

//partial specialization for double
template<>
struct PythonUtility::convertFromPython<double>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static double get(PyObject *object, double defaultValue)
  {
    if(object == NULL)
      return defaultValue;

    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    if (PyFloat_Check(object))
    {
      double valueDouble = PyFloat_AsDouble(object);

      return valueDouble;
    }
    else if (PyLong_Check(object))
    {
      long valueLong = PyLong_AsLong(object);
      return double(valueLong);
    }
    else if (PyUnicode_Check(object))
    {
      std::string valueString = pyUnicodeToString(object);
      return atof(valueString.c_str());
    }
    else if (PyComplex_Check(object))
    {
      return PyComplex_RealAsDouble(object);
    }
    /*
  #ifdef HAVE_NUMPYC
    else if (PyArray_Check(object))
    {
      //if (object->descr->type_num != NPY_DOUBLE || vec->nd != 1)  {

      LOG(WARNING) << "convertFromPython: object is a numpy array: " << object;
    }
  #endif
  */
    else
    {
      LOG(WARNING) << "convertFromPython: object is no double: " << object;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static double get(PyObject *object)
  {
    return convertFromPython<double>::get(object, 0.0);
  }
};

//partial specialization for std::string
template<>
struct PythonUtility::convertFromPython<std::string>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static std::string get(PyObject *object, std::string defaultValue)
  {
    if(object == NULL)
      return defaultValue;

    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    if (PyUnicode_Check(object))
    {
      std::string valueString = pyUnicodeToString(object);
      return valueString;
    }
    else
    {
      LOG(WARNING) << "convertFromPython: object is no std::string: " << object;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static std::string get(PyObject *object)
  {
    return convertFromPython<std::string>::get(object, "");
  }
};

//partial specialization for global_no_t
template<>
struct PythonUtility::convertFromPython<global_no_t>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static global_no_t get(PyObject *object, global_no_t defaultValue)
  {
    if(object == NULL)
      return defaultValue;

    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    if (PyLong_Check(object))
    {
      long valueLong = PyLong_AsLong(object);
      return global_no_t(valueLong);
    }
    else if (PyFloat_Check(object))
    {
      double valueDouble = PyFloat_AsDouble(object);

      if (double(global_no_t(valueDouble)) != valueDouble)      // if value is not e.g. 2.0
      {
        LOG(WARNING) << "convertFromPython: object is no int: " << object;
      }

      return global_no_t(valueDouble);
    }
    else if (PyUnicode_Check(object))
    {
      std::string valueString = pyUnicodeToString(object);
      return atoi(valueString.c_str());
    }
    else
    {
      LOG(WARNING) << "convertFromPython: object is no long int: " << object;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static global_no_t get(PyObject *object)
  {
    return convertFromPython<global_no_t>::get(object, 0);
  }
};

//partial specialization for std::size_t
template<>
struct PythonUtility::convertFromPython<std::size_t>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static std::size_t get(PyObject *object, std::size_t defaultValue)
  {
    if(object == NULL)
      return defaultValue;

    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    if (PyLong_Check(object))
    {
      long valueLong = PyLong_AsLong(object);
      return std::size_t(valueLong);
    }
    else if (PyFloat_Check(object))
    {
      double valueDouble = PyFloat_AsDouble(object);

      if (double(std::size_t(valueDouble)) != valueDouble)      // if value is not e.g. 2.0
      {
        LOG(WARNING) << "convertFromPython: object is no std::size_t: " << object;
      }

      return std::size_t(valueDouble);
    }
    else if (PyUnicode_Check(object))
    {
      std::string valueString = pyUnicodeToString(object);
      return atoi(valueString.c_str());
    }
    else
    {
      LOG(WARNING) << "convertFromPython: object is no std::size_t: " << object;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static std::size_t get(PyObject *object)
  {
    return convertFromPython<std::size_t>::get(object, 0);
  }
};

//partial specialization for PyObject *
template<>
struct PythonUtility::convertFromPython<PyObject *>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static PyObject *get(PyObject *object, PyObject *defaultValue)
  {
    return object;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static PyObject *get(PyObject *object)
  {
    return convertFromPython<PyObject *>::get(object, nullptr);
  }

};

//partial specialization for bool
template<>
struct PythonUtility::convertFromPython<bool>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static bool get(PyObject *object, bool defaultValue)
  {
    if(object == NULL)
      return defaultValue;

    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    if (PyBool_Check(object))
    {
      if (object == Py_True)
      {
        return true;
      }
      return false;
    }
    if (PyLong_Check(object))
    {
      long valueLong = PyLong_AsLong(object);
      return bool(valueLong);
    }
    else if (PyUnicode_Check(object))
    {
      std::string valueString = pyUnicodeToString(object);
      std::transform(valueString.begin(), valueString.end(), valueString.begin(), ::tolower);
      if (valueString.find("true") != std::string::npos || valueString.find("1") != std::string::npos
        || valueString.find("yes") != std::string::npos || valueString.find("on") != std::string::npos)
      {
        return true;
      }
      else if (valueString.find("false") != std::string::npos || valueString.find("0") != std::string::npos
        || valueString.find("no") != std::string::npos || valueString.find("off") != std::string::npos)
      {
        return false;
      }
      else
      {
        LOG(WARNING) << "Could not infer bool value of \"" <<valueString<< "\".";
        return defaultValue;
      }
    }
    else
    {
      LOG(WARNING) << "convertFromPython: object is no bool: " << object;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static bool get(PyObject *object)
  {
    return convertFromPython<bool>::get(object, false);
  }
};

// ------------- convertToPython ----------------
template<>
struct PythonUtility::convertToPython<int>
{
  //! convert a type to a python object
  static PyObject *get(int value)
  {
    return PyLong_FromLong(value);
  }
};

template<>
struct PythonUtility::convertToPython<double>
{
  //! convert a type to a python object
  static PyObject *get(double value)
  {
    return PyFloat_FromDouble(value);
  }
};

template<>
struct PythonUtility::convertToPython<std::string>
{
  //! convert a type to a python object
  static PyObject *get(std::string value)
  {
    return PyUnicode_FromString(value.c_str());
  }
};

template<>
struct PythonUtility::convertToPython<global_no_t>
{
  //! convert a type to a python object
  static PyObject *get(global_no_t value)
  {
    return PyLong_FromLongLong(value);
  }
};

template<>
struct PythonUtility::convertToPython<std::size_t>
{
  //! convert a type to a python object
  static PyObject *get(std::size_t value)
  {
    return PyLong_FromLongLong(value);
  }
};

template<>
struct PythonUtility::convertToPython<PyObject *>
{
  //! convert a type to a python object
  static PyObject *get(PyObject * value)
  {
    return value;
  }
};

template<>
struct PythonUtility::convertToPython<bool>
{
  //! convert a type to a python object
  static PyObject *get(bool value)
  {
    return PyBool_FromLong(value);
  }
};

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
