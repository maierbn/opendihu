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

//partial specialization for int
template<>
struct PythonUtility::convertFromPython<long>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static int get(PyObject *object, long defaultValue)
  {
    if (object == NULL)
      return defaultValue;

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    assert(object != nullptr);
    if (PyLong_Check(object))
    {
      long valueLong = PyLong_AsLong(object);
      return valueLong;
    }
    else if (PyFloat_Check(object))
    {
      double valueDouble = PyFloat_AsDouble(object);

      if (double(int(valueDouble)) != valueDouble)      // if value is not e.g. 2.0
      {
        LOG(WARNING) << "convertFromPython<long>: object is float and not long: " << object;
      }

      return int(valueDouble);
    }
    else if (PyUnicode_Check(object))
    {
      std::string valueString = pyUnicodeToString(object);
      return atoi(valueString.c_str());
    }
    else if (object == Py_None)
    {
      LOG(DEBUG) << "convertFromPython<long>: object is None, parse as -1";
      return -1;    // None translates to -1
    }
    else
    {
      LOG(WARNING) << "convertFromPython<long>: object is no long: " << object;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static int get(PyObject *object)
  {
    return convertFromPython<long>::get(object, 0);
  }
};

//partial specialization for int
template<>
struct PythonUtility::convertFromPython<int>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static int get(PyObject *object, int defaultValue)
  {
    if (object == NULL)
      return defaultValue;

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    assert(object != nullptr);
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
        LOG(WARNING) << "convertFromPython<int>: object is float and not int: " << object;
      }

      return int(valueDouble);
    }
    else if (PyUnicode_Check(object))
    {
      std::string valueString = pyUnicodeToString(object);
      return atoi(valueString.c_str());
    }
    else if (object == Py_None)
    {
      LOG(DEBUG) << "convertFromPython<int>: object is None, parse as -1";
      return -1;    // None translates to -1
    }
    else
    {
      LOG(WARNING) << "convertFromPython<int>: object is no int: " << object;
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
    if (object == NULL)
      return defaultValue;

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    if (PyFloat_Check(object))
    {
      double valueDouble = PyFloat_AsDouble(object);

      return valueDouble;
    }
    else if (object == Py_None)
    {
      //return std::nan("");
      return std::numeric_limits<double>::max();
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
    if (object == NULL)
      return defaultValue;

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    if (PyUnicode_Check(object))
    {
      std::string valueString = pyUnicodeToString(object);
      return valueString;
    }
    else
    {
      PyObject *unicode = PyObject_Str(object);

      if (unicode != NULL)
      {
        std::string valueString = pyUnicodeToString(unicode);
        return valueString;
      }

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
    if (object == NULL)
      return defaultValue;

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

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
        LOG(WARNING) << "convertFromPython: object is float and not int or global_no_t: " << object;
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
    if (object == NULL)
      return defaultValue;

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

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

//partial specialization for long long which is also MPI_Offset
template<>
struct PythonUtility::convertFromPython<long long>
{
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  static long long get(PyObject *object, long long defaultValue)
  {
    if (object == NULL)
      return defaultValue;

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    if (PyLong_Check(object))
    {
      long long valueLong = PyLong_AsLongLong(object);
      return valueLong;
    }
    else if (PyFloat_Check(object))
    {
      double valueDouble = PyFloat_AsDouble(object);

      if (double((long long)(valueDouble)) != valueDouble)      // if value is not e.g. 2.0
      {
        LOG(WARNING) << "convertFromPython: object is no long long: " << object;
      }

      return (long long)(valueDouble);
    }
    else if (PyUnicode_Check(object))
    {
      std::string valueString = pyUnicodeToString(object);
      return atoi(valueString.c_str());
    }
    else
    {
      LOG(WARNING) << "convertFromPython: object is no long long: " << object;
    }
    return defaultValue;
  }

  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  static long long get(PyObject *object)
  {
    return convertFromPython<long long>::get(object, 0);
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
    if (object == NULL)
      return defaultValue;

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

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
