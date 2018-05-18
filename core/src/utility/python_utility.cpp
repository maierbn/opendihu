#include "utility/python_utility.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <Python.h>
#include "easylogging++.h"

#include "control/use_numpy.h"

PyObject *PythonUtility::itemList = NULL;
int PythonUtility::itemListIndex = 0;
PyObject *PythonUtility::list = NULL;
int PythonUtility::listIndex = 0;

template<>
int PythonUtility::convertFromPython(PyObject *object, int defaultValue)
{
  if(object == NULL)
    return defaultValue;

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

template<>
std::size_t PythonUtility::convertFromPython(PyObject *object, std::size_t defaultValue)
{
  if(object == NULL)
    return defaultValue;

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

template<>
double PythonUtility::convertFromPython(PyObject *object, double defaultValue)
{
  initNumpy();

  if(object == NULL)
    return defaultValue;

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

template<>
std::string PythonUtility::convertFromPython(PyObject *object, std::string defaultValue)
{
  if(object == NULL)
    return defaultValue;

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

template<>
PyObject *PythonUtility::convertFromPython(PyObject *object, PyObject *defaultValue)
{
  return object;
}

template<>
int PythonUtility::convertFromPython(PyObject *object)
{
  return convertFromPython<int>(object, 0);
}

template<>
std::size_t PythonUtility::convertFromPython(PyObject *object)
{
  return convertFromPython<std::size_t>(object, 0);
}

template<>
bool PythonUtility::convertFromPython(PyObject *object, bool defaultValue)
{
  if(object == NULL)
    return defaultValue;

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
      LOG(WARNING) << "Could not infer bool value of \""<<valueString<<"\".";
      return defaultValue;
    }
  }
  else
  {
    LOG(WARNING) << "convertFromPython: object is no bool: " << object;
  }
  return defaultValue;
}

template<>
double PythonUtility::convertFromPython(PyObject *object)
{
  return convertFromPython<double>(object, 0.0);
}

template<>
std::string PythonUtility::convertFromPython(PyObject *object)
{
  return convertFromPython<std::string>(object, "");
}

template<>
PyObject *PythonUtility::convertFromPython(PyObject *object)
{
  return convertFromPython<PyObject *>(object, NULL);
}

template<>
std::array<double,2> PythonUtility::convertFromPython(PyObject *object, std::array<double,2> defaultValue)
{
  return PythonUtility::convertFromPython<double,2>(object, defaultValue);
}

template<>
std::array<double,2> PythonUtility::convertFromPython(PyObject *object)
{
  return convertFromPython<double,2>(object);
}

template<>
std::array<double,3> PythonUtility::convertFromPython(PyObject *object, std::array<double,3> defaultValue)
{
  return PythonUtility::convertFromPython<double,3>(object, defaultValue);
}

template<>
std::array<double,3> PythonUtility::convertFromPython(PyObject *object)
{
  return convertFromPython<double,3>(object);
}

bool PythonUtility::hasKey(const PyObject* settings, std::string keyString)
{
  if (settings)
  {
    // check if input dictionary contains the key
    PyObject *key = PyUnicode_FromString(keyString.c_str());

    if(PyDict_Contains((PyObject *)settings, key))
    {
      Py_CLEAR(key);
      return true;
    }
    Py_CLEAR(key);
  }
  return false;
}

PyObject *PythonUtility::getOptionPyObject(const PyObject *settings, std::string keyString)
{
  if (settings)
  {
    // check if input dictionary contains the key
    PyObject *key = PyUnicode_FromString(keyString.c_str());
    if(PyDict_Contains((PyObject *)settings, key))
    {
      PyObject *value = PyDict_GetItem((PyObject *)settings, key);
      Py_CLEAR(key);
      return value;
    }
    else
    {
      LOG(WARNING) << "Dict does not contain Key \""<<keyString<<"\"!"<<std::endl;
      Py_CLEAR(key);
      return NULL;
    }
  }
  return NULL;
}

double PythonUtility::getOptionDouble(const PyObject* settings, std::string keyString, double defaultValue, ValidityCriterion validityCriterion)
{
  double result = defaultValue;

  if (!settings)
  {
    LOG(DEBUG)<<"PyObject *settings is NULL.";
    return result;
  }

  // check if input dictionary contains the key
  PyObject *key = PyUnicode_FromString(keyString.c_str());
  if(PyDict_Contains((PyObject *)settings, key))
  {
    // extract the value of the key and check its type
    PyObject *value = PyDict_GetItem((PyObject *)settings, key);
    if (PyList_Check(value))
    {
      // type is a list, determine how many entries it contains
      int listNEntries = PyList_Size(value);
      LOG(DEBUG)<<"list with "<<listNEntries<<" entries"<<std::endl;

      // if there are multiple entries, use the first
      if (listNEntries >= 1)
      {
        // extract first entry
        PyObject *listEntry = PyList_GetItem(value, (Py_ssize_t)0);

        // convert to double
        result = convertFromPython<double>(listEntry, defaultValue);

        // print a warning if there are further entries
        if (listNEntries > 1)
        {
          LOG(WARNING)<<"Only using first value "<<result<<" of list of length "<<listNEntries
                        <<" at Key \""<<keyString<<"\"."<<std::endl;
        }
      }
      else
      {
        LOG(WARNING)<<"Empty list for Key \""<<keyString<<"\"."<<std::endl;
      }
    }
    else
    {
      // convert to double or take default value
      result = convertFromPython<double>(value, defaultValue);

      //LOG(DEBUG)<<"PythonUtility::getOptionDouble: Value for key \""<<keyString<<"\" found: "<<result<<".";
    }
  }
  else
  {
    LOG(WARNING)<<"Key \""<<keyString<<"\" not found in dict in config file.";
  }

  switch(validityCriterion)
  {
  case Positive:
    if (result <= 0.0)
    {
      LOG(WARNING)<<"value "<<result<<" of Key \""<<keyString<<"\" is invalid (not positive). Using default value "
        <<defaultValue<<".";
      result = defaultValue;
    }
    break;
  case NonNegative:
    if (result < 0.0)
    {
      LOG(WARNING)<<"value "<<result<<" of Key \""<<keyString<<"\" is invalid (not non-negative). Using default value "
        <<defaultValue<<".";
      result = defaultValue;
    }
    break;
  case Between1And3:
    if (result < 1.0)
    {
      LOG(WARNING)<<"value "<<result<<" of Key \""<<keyString<<"\" is invalid (<1). Using default value "
        <<defaultValue<<".";
      result = defaultValue;
    }
    else if (result > 3.0)
    {
      LOG(WARNING)<<"value "<<result<<" of Key \""<<keyString<<"\" is invalid (>3). Using default value "
        <<defaultValue<<".";
      result = defaultValue;
    }
    break;
  case None:
    break;
  };

  Py_CLEAR(key);
  return result;
}

int PythonUtility::getOptionInt(const PyObject *settings, std::string keyString, int defaultValue, ValidityCriterion validityCriterion)
{
  int result = defaultValue;

  if (!settings)
    return result;

  // check if input dictionary contains the key
  PyObject *key = PyUnicode_FromString(keyString.c_str());
  if(PyDict_Contains((PyObject *)settings, key))
  {
    // extract the value of the key and check its type
    PyObject *value = PyDict_GetItem((PyObject *)settings, key);
    if (PyList_Check(value))
    {
      // type is a list, determine how many entries it contains
      int listNEntries = PyList_Size(value);
      LOG(DEBUG)<<"list with "<<listNEntries<<" entries"<<std::endl;

      // if there are multiple entries, use the first
      if (listNEntries >= 1)
      {
        // extract first entry
        PyObject *listEntry = PyList_GetItem(value, (Py_ssize_t)0);

        // convert to int
        result = convertFromPython<int>(listEntry, defaultValue);

        // print a warning if there are further entries
        if (listNEntries > 1)
        {
          LOG(WARNING)<<"Only using first value "<<result<<" of list of length "<<listNEntries
            <<" at Key \""<<keyString<<"\"."<<std::endl;
        }
      }
      else
      {
        LOG(WARNING)<<"Empty list for Key \""<<keyString<<"\"."<<std::endl;
      }
    }
    else
    {
      // convert to int or try default value
      result = convertFromPython<int>(value, defaultValue);
    }
  }
  else
  {
    LOG(WARNING)<<"Key \""<<keyString<<"\" not found in dict in config file.";
  }

  switch(validityCriterion)
  {
    case Positive:
      if (result <= 0)
      {
        LOG(WARNING)<<"value "<<result<<" of Key \""<<keyString<<"\" is invalid (not positive). Using default value "
          <<defaultValue<<".";
        result = defaultValue;
      }
      break;
    case NonNegative:
      if (result < 0)
      {
        LOG(WARNING)<<"value "<<result<<" of Key \""<<keyString<<"\" is invalid (not non-negative). Using default value "
          <<defaultValue<<".";
        result = defaultValue;
      }
      break;
    case Between1And3:
      if (result < 1)
      {
        LOG(WARNING)<<"value "<<result<<" of Key \""<<keyString<<"\" is invalid (<1). Using default value "
          <<defaultValue<<".";
        result = defaultValue;
      }
      else if (result > 3)
      {
        LOG(WARNING)<<"value "<<result<<" of Key \""<<keyString<<"\" is invalid (>3). Using default value "
          <<defaultValue<<".";
        result = defaultValue;
      }
      break;
    case None:
      break;
  };

  Py_CLEAR(key);
  return result;
}

bool PythonUtility::getOptionBool(const PyObject *settings, std::string keyString, bool defaultValue)
{
  int result = defaultValue;

  if (!settings)
    return result;

  // check if input dictionary contains the key
  PyObject *key = PyUnicode_FromString(keyString.c_str());
  if(PyDict_Contains((PyObject *)settings, key))
  {
    // extract the value of the key and check its type
    PyObject *value = PyDict_GetItem((PyObject *)settings, key);
    if (PyList_Check(value))
    {
      // type is a list, determine how many entries it contains
      int listNEntries = PyList_Size(value);
      LOG(DEBUG)<<"list with "<<listNEntries<<" entries"<<std::endl;

      // if there are multiple entries, use the first
      if (listNEntries >= 1)
      {
        // extract first entry
        PyObject *listEntry = PyList_GetItem(value, (Py_ssize_t)0);

        // convert to bool
        result = convertFromPython<bool>(listEntry, defaultValue);

        // print a warning if there are further entries
        if (listNEntries > 1)
        {
          LOG(WARNING)<<"Only using first value "<<result<<" of list of length "<<listNEntries
            <<" at Key \""<<keyString<<"\"."<<std::endl;
        }
      }
      else
      {
        LOG(WARNING)<<"Empty list for Key \""<<keyString<<"\"."<<std::endl;
      }
    }
    else
    {
      // convert to bool or try default value
      result = convertFromPython<bool>(value, defaultValue);
    }
  }
  else
  {
    LOG(WARNING)<<"Key \""<<keyString<<"\" not found in dict in config file.";
  }
  Py_CLEAR(key);
  return result;
}

std::string PythonUtility::getOptionString(const PyObject *settings, std::string keyString, std::string defaultValue)
{
  std::string result = defaultValue;

  if (!settings)
    return result;

  // check if input dictionary contains the key
  PyObject *key = PyUnicode_FromString(keyString.c_str());
  if(PyDict_Contains((PyObject *)settings, key))
  {
    // extract the value of the key and check its type
    PyObject *value = PyDict_GetItem((PyObject *)settings, key);

    // convert to std::string or try default value
    result = convertFromPython<std::string>(value, defaultValue);
  }
  else
  {
    LOG(WARNING)<<"Key \""<<keyString<<"\" not found in dict in config file.";
  }

  Py_CLEAR(key);
  return result;
}

PyObject *PythonUtility::getOptionFunction(const PyObject *settings, std::string keyString)
{
  PyObject *result = NULL;

  if (!settings)
    return result;

  // check if input dictionary contains the key
  PyObject *key = PyUnicode_FromString(keyString.c_str());
  if(PyDict_Contains((PyObject *)settings, key))
  {
    // extract the value of the key and check its type
    PyObject *function = PyDict_GetItem((PyObject *)settings, key);
    if (PyFunction_Check(function))
    {
      // type is a function
      if (PyCallable_Check(function))
      {
        result = function;
      }
      else
      {
        LOG(WARNING)<<"Value for key \""<<keyString<<"\" is not a callable object.";
      }
    }
    else
    {
      LOG(WARNING)<<"Value for key \""<<keyString<<"\" is not a function.";
    }
  }
  else
  {
    LOG(WARNING)<<"Key \""<<keyString<<"\" not found in dict in config file.";
  }
  Py_CLEAR(key);
  return result;
}

std::string PythonUtility::getString(PyObject *object, int indent, int first_indent)
{
  if (!object)
    return "NULL";

  std::stringstream line;
  line << std::string(first_indent, ' ');

  if (PyUnicode_CheckExact(object))
  {
    std::string objectString = pyUnicodeToString(object);
    line << "\""<<objectString<<"\"";
  }
  else if (PyLong_CheckExact(object))
  {
    long objectLong = PyLong_AsLong(object);
    line << objectLong;
  }
  else if (PyLong_CheckExact(object))
  {
    long objectLong = PyLong_AsLong(object);
    line << objectLong;
  }
  else if (PyFloat_CheckExact(object))
  {
    double objectDouble = PyFloat_AsDouble(object);
    line << objectDouble;
  }
  else if (PyBool_Check(object))
  {
    bool objectBool = PyObject_IsTrue(object);
    line << std::boolalpha<<objectBool;
  }
  else if (PyList_Check(object))
  {
    int size = (int)PyList_Size(object);
    line << "[";
    for (int index = 0; index < size; index++)
    {
      PyObject *item = PyList_GetItem(object, (Py_ssize_t)index);

      line << getString(item, indent+2, 0) << (index < size-1? ", " : "]");
    }
  }
  else if(PyDict_CheckExact(object))
  {
    // iterate over top level key-value pairs
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    line << "{";
    bool first = true;
    while (PyDict_Next(object, &pos, &key, &value))
    {
      if (!first)
        line << ",";
      first = false;

      line << std::endl << std::string(indent+2, ' ');

      if (PyUnicode_Check(key))
      {
        std::string keyString = pyUnicodeToString(key);
        line << keyString<<": ";
      }
      else if (PyLong_Check(key))
      {
        std::string keyString = std::to_string(PyLong_AsLong(key));
        line << keyString<<": ";
      }
      else
      {
        line << "(key is of unknown type): ";
      }

      line << getString(value, indent+2, 0);
    }
    line << std::endl << std::string(indent, ' ') << "}";
  }
  else
  {
    line << "<unknown type>";
  }

  return line.str();
}

void PythonUtility::printDict(PyObject *dict)
{
  if (dict == NULL)
  {
    VLOG(1) << "dict is NULL!";
    return;
  }

  if (!PyDict_Check(dict))
  {
    VLOG(1) << "Object is not a dict!";
    return;
  }

  VLOG(1) << getString(dict);
}

bool PythonUtility::getOptionDictEnd(const PyObject *settings, std::string keyString)
{
  if (!itemList)
    return true;
  return itemListIndex >= PyList_Size(itemList);
}

bool PythonUtility::getOptionListEnd(const PyObject *settings, std::string keyString)
{
  if (!list)
    return true;
  return listIndex >= PyList_Size(list);
}

void PythonUtility::getOptionVector(const PyObject* settings, std::string keyString, int nEntries, std::vector<double> &values)
{
  values.resize(nEntries);

  if (settings)
  {

    // check if input dictionary contains the key
    PyObject *key = PyUnicode_FromString(keyString.c_str());
    if(PyDict_Contains((PyObject *)settings, key))
    {
      // extract the value of the key and check its type
      PyObject *value = PyDict_GetItem((PyObject *)settings, key);
      if (PyList_Check(value))
      {
        // it is a list

        // get the first value from the list
        double value = PythonUtility::getOptionListBegin<double>(settings, keyString);
        int i = 0;

        // loop over other values
        for (;
            !PythonUtility::getOptionListEnd(settings, keyString)
            && i < nEntries;
            PythonUtility::getOptionListNext<double>(settings, keyString, value), i++)
        {
          values[i] = value;
        }
      }
      else if (PyDict_Check(value))
      {
        std::pair<int, double> dictItem
          = PythonUtility::getOptionDictBegin<int, double>(settings, keyString);

        // loop over Dirichlet boundary conditions
        for (; !PythonUtility::getOptionDictEnd(settings, keyString);
            PythonUtility::getOptionDictNext<int, double>(settings, keyString, dictItem))
        {
          int index = dictItem.first;
          double value = dictItem.second;

          if (index >= 0 && index < nEntries)
          {
            values[index] = value;
          }
          else
          {
            LOG(WARNING) << "In config dict, ignoring key "<<index<<", maximum key is "<<nEntries-1;
          }
        }
      }
      else
      {
        double value = PythonUtility::getOptionDouble(settings, keyString, 0.0);
        std::fill(values.begin(), values.end(), value);
      }
    }
    else
    {
      LOG(WARNING)<<"Key \""<<keyString<<"\" not found in dict in config file.";
    }
    Py_CLEAR(key);
  }
}

PyObject *PythonUtility::convertToPythonList(std::vector<double> &data)
{
  PyObject *result = PyList_New((Py_ssize_t)data.size());
  for (unsigned int i=0; i<data.size(); i++)
  {
    PyObject *item = PyFloat_FromDouble(data[i]);
    PyList_SetItem(result, (Py_ssize_t)i, item);    // steals reference to item
  }
  return result;    // return value: new reference
}

PyObject *PythonUtility::convertToPythonList(std::vector<long> &data)
{
  PyObject *result = PyList_New((Py_ssize_t)data.size());
  for (unsigned int i=0; i<data.size(); i++)
  {
    PyObject *item = PyLong_FromLong(data[i]);
    PyList_SetItem(result, (Py_ssize_t)i, item);    // steals reference to item
  }
  return result;    // return value: new reference
}

PyObject *PythonUtility::convertToPythonList(unsigned int nEntries, double* data)
{
  PyObject *result = PyList_New((Py_ssize_t)nEntries);
  for (unsigned int i=0; i<nEntries; i++)
  {
    PyObject *item = PyFloat_FromDouble(data[i]);
    PyList_SetItem(result, (Py_ssize_t)i, item);    // steals reference to item
  }
  return result;    // return value: new reference
}

std::string PythonUtility::pyUnicodeToString(PyObject* object)
{
#if PY_MAJOR_VERSION >= 3
  PyObject *asciiString = PyUnicode_AsASCIIString(object);
  std::string result = PyBytes_AsString(asciiString);
  Py_DECREF(asciiString);
#else
  std::string result = pyUnicodeToString(object);
#endif

  return result;
}

std::ostream &operator<<(std::ostream &stream, PyObject *object)
{
  stream << PythonUtility::getString(object);
  return stream;
}