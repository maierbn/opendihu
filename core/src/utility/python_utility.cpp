#include "utility/python_utility.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <mutex>          // std::mutex, std::unique_lock, std::defer_lock
#include <omp.h>

#include <Python.h>
#include "easylogging++.h"

#include "control/use_numpy.h"
#include "control/types.h"

PyObject *PythonUtility::itemList = NULL;
int PythonUtility::itemListIndex = 0;
PyObject *PythonUtility::list = NULL;
int PythonUtility::listIndex = 0;

bool PythonUtility::hasKey(const PyObject* settings, std::string keyString)
{
  if (settings)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;
  
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

bool PythonUtility::isTypeList(const PyObject *object)
{
  if (object)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;
  
    if (PyList_Check(object))
    {
      return true;
    }
  }
  return false;
}

PyObject *PythonUtility::getOptionPyObject(const PyObject *settings, std::string keyString)
{
  if (settings)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;
    
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
      LOG(WARNING) << "Dict does not contain Key \"" <<keyString<< "\"!" << std::endl;
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
    LOG(DEBUG) << "PyObject *settings is NULL.";
    return result;
  }

  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
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
      LOG(DEBUG) << "list with " << listNEntries << " entries" << std::endl;

      // if there are multiple entries, use the first
      if (listNEntries >= 1)
      {
        // extract first entry
        PyObject *listEntry = PyList_GetItem(value, (Py_ssize_t)0);

        // convert to double
        result = convertFromPython<double>::get(listEntry, defaultValue);

        // print a warning if there are further entries
        if (listNEntries > 1)
        {
          LOG(WARNING) << "Only using first value " <<result<< " of list of length " << listNEntries
                        << " at Key \"" <<keyString<< "\"." << std::endl;
        }
      }
      else
      {
        LOG(WARNING) << "Empty list for Key \"" <<keyString<< "\"." << std::endl;
      }
    }
    else
    {
      // convert to double or take default value
      result = convertFromPython<double>::get(value, defaultValue);

      //LOG(DEBUG) << "PythonUtility::getOptionDouble: Value for key \"" <<keyString<< "\" found: " <<result<< ".";
    }
  }
  else
  {
    LOG(WARNING) << "Key \"" <<keyString<< "\" not found in dict in config file. Assuming default value " << defaultValue << ".";
  }

  switch(validityCriterion)
  {
  case Positive:
    if (result <= 0.0)
    {
      LOG(WARNING) << "value " <<result<< " of Key \"" <<keyString<< "\" is invalid (not positive). Using default value "
        << defaultValue << ".";
      result = defaultValue;
    }
    break;
  case NonNegative:
    if (result < 0.0)
    {
      LOG(WARNING) << "value " <<result<< " of Key \"" <<keyString<< "\" is invalid (not non-negative). Using default value "
        << defaultValue << ".";
      result = defaultValue;
    }
    break;
  case Between1And3:
    if (result < 1.0)
    {
      LOG(WARNING) << "value " <<result<< " of Key \"" <<keyString<< "\" is invalid (<1). Using default value "
        << defaultValue << ".";
      result = defaultValue;
    }
    else if (result > 3.0)
    {
      LOG(WARNING) << "value " <<result<< " of Key \"" <<keyString<< "\" is invalid (>3). Using default value "
        << defaultValue << ".";
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

  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
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
      LOG(DEBUG) << "list with " << listNEntries << " entries" << std::endl;

      // if there are multiple entries, use the first
      if (listNEntries >= 1)
      {
        // extract first entry
        PyObject *listEntry = PyList_GetItem(value, (Py_ssize_t)0);

        // convert to int
        result = convertFromPython<int>::get(listEntry, defaultValue);

        // print a warning if there are further entries
        if (listNEntries > 1)
        {
          LOG(WARNING) << "Only using first value " <<result<< " of list of length " << listNEntries
            << " at Key \"" <<keyString<< "\"." << std::endl;
        }
      }
      else
      {
        LOG(WARNING) << "Empty list for Key \"" <<keyString<< "\"." << std::endl;
      }
    }
    else
    {
      // convert to int or try default value
      result = convertFromPython<int>::get(value, defaultValue);
    }
  }
  else
  {
    LOG(WARNING) << "Key \"" <<keyString<< "\" not found in dict in config file. Assuming default value " << defaultValue << ".";
  }

  switch(validityCriterion)
  {
    case Positive:
      if (result <= 0)
      {
        LOG(WARNING) << "value " <<result<< " of Key \"" <<keyString<< "\" is invalid (not positive). Using default value "
          << defaultValue << ".";
        result = defaultValue;
      }
      break;
    case NonNegative:
      if (result < 0)
      {
        LOG(WARNING) << "value " <<result<< " of Key \"" <<keyString<< "\" is invalid (not non-negative). Using default value "
          << defaultValue << ".";
        result = defaultValue;
      }
      break;
    case Between1And3:
      if (result < 1)
      {
        LOG(WARNING) << "value " <<result<< " of Key \"" <<keyString<< "\" is invalid (<1). Using default value "
          << defaultValue << ".";
        result = defaultValue;
      }
      else if (result > 3)
      {
        LOG(WARNING) << "value " <<result<< " of Key \"" <<keyString<< "\" is invalid (>3). Using default value "
          << defaultValue << ".";
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

  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
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
      LOG(DEBUG) << "list with " << listNEntries << " entries" << std::endl;

      // if there are multiple entries, use the first
      if (listNEntries >= 1)
      {
        // extract first entry
        PyObject *listEntry = PyList_GetItem(value, (Py_ssize_t)0);

        // convert to bool
        result = convertFromPython<bool>::get(listEntry, defaultValue);

        // print a warning if there are further entries
        if (listNEntries > 1)
        {
          LOG(WARNING) << "Only using first value " <<result<< " of list of length " << listNEntries
            << " at Key \"" <<keyString<< "\"." << std::endl;
        }
      }
      else
      {
        LOG(WARNING) << "Empty list for Key \"" <<keyString<< "\"." << std::endl;
      }
    }
    else
    {
      // convert to bool or try default value
      result = convertFromPython<bool>::get(value, defaultValue);
    }
  }
  else
  {
    LOG(WARNING) << "Key \"" <<keyString<< "\" not found in dict in config file. Assuming default value " << std::boolalpha << defaultValue << ".";
  }
  Py_CLEAR(key);
  return result;
}

std::string PythonUtility::getOptionString(const PyObject *settings, std::string keyString, std::string defaultValue)
{
  std::string result = defaultValue;

  if (!settings)
    return result;

  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
  // check if input dictionary contains the key
  PyObject *key = PyUnicode_FromString(keyString.c_str());
  if(PyDict_Contains((PyObject *)settings, key))
  {
    // extract the value of the key and check its type
    PyObject *value = PyDict_GetItem((PyObject *)settings, key);

    // convert to std::string or try default value
    result = convertFromPython<std::string>::get(value, defaultValue);
  }
  else
  {
    LOG(WARNING) << "Key \"" <<keyString<< "\" not found in dict in config file. Assuming default value \"" << defaultValue << "\".";
  }

  Py_CLEAR(key);
  return result;
}

PyObject *PythonUtility::getOptionFunction(const PyObject *settings, std::string keyString)
{
  PyObject *result = NULL;

  if (!settings)
    return result;

  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
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
        LOG(WARNING) << "Value for key \"" <<keyString<< "\" is not a callable object.";
      }
    }
    else
    {
      LOG(WARNING) << "Value for key \"" <<keyString<< "\" is not a function.";
    }
  }
  else
  {
    LOG(WARNING) << "Key \"" <<keyString<< "\" not found in dict in config file.";
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

  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
  if (PyUnicode_CheckExact(object))
  {
    std::string objectString = pyUnicodeToString(object);
    line << "\"" <<objectString<< "\"";
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
  else if (PyComplex_Check(object))
  {
    double realDouble = PyComplex_RealAsDouble(object);
    double imagDouble = PyComplex_ImagAsDouble(object);
    line << realDouble << " + " << imagDouble << "i";
  }
  else if (PyBytes_Check(object))
  {
    std::string str = PyBytes_AsString(object);
    line << "\"" << str << "\"";
  }
  else if (PyBool_Check(object))
  {
    bool objectBool = PyObject_IsTrue(object);
    line << std::boolalpha << objectBool;
  }
  else if (PyUnicode_Check(object))
  {
    std::string str = PyUnicode_AS_DATA(object);
    line << "\"" << str << "\"";
  }
  else if (PyFunction_Check(object))
  {
    line << "<function object>";
  }
  else if (object == Py_True)
  {
    line << "True";
  }
  else if (object == Py_False)
  {
    line << "False";
  }
  else if (object == Py_None)
  {
    line << "None";
  }
  else if (PyCallIter_Check(object))
  {
    line << "<iterator object>";
  }
  else if (PyByteArray_Check(object))
  {
    line << "<byte array object>";
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
  else if (PyTuple_Check(object))
  {
    int size = (int)PyTuple_Size(object);
    line << "(";
    for (int index = 0; index < size; index++)
    {
      PyObject *item = PyTuple_GetItem(object, (Py_ssize_t)index);

      line << getString(item, indent+2, 0) << (index < size-1? ", " : ")");
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
        line << keyString<< ": ";
      }
      else if (PyLong_Check(key))
      {
        std::string keyString = std::to_string(PyLong_AsLong(key));
        line << keyString<< ": ";
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
    LOG(DEBUG) << "printDict: dict is NULL!";
    return;
  }

  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
  if (!PyDict_Check(dict))
  {
    LOG(DEBUG) << "printDict: Object is not a dict!";
    return;
  }

  LOG(DEBUG) << getString(dict);
}

bool PythonUtility::getOptionDictEnd(const PyObject *settings, std::string keyString)
{
  if (!itemList)
    return true;
  
  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
  return itemListIndex >= PyList_Size(itemList);
}

bool PythonUtility::getOptionListEnd(const PyObject *settings, std::string keyString)
{
  if (!list)
    return true;
  
  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
  return listIndex >= PyList_Size(list);
}

void PythonUtility::getOptionVector(const PyObject* settings, std::string keyString, int nEntries, std::vector<double> &values)
{
  values.resize(nEntries, 0.0);

  if (settings)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;
  
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
            LOG(WARNING) << "In config dict, ignoring key " <<index<< ", maximum key is " <<nEntries-1;
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
      LOG(WARNING) << "Key \"" <<keyString<< "\" not found in dict in config file, assuming list with " << nEntries << " zeros";
    }
    Py_CLEAR(key);
  }
}

void PythonUtility::getOptionVector(const PyObject *settings, std::string keyString, std::vector<int> &values)
{
  if (settings)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    // check if input dictionary contains the key
    PyObject *key = PyUnicode_FromString(keyString.c_str());
    if(PyDict_Contains((PyObject *)settings, key))
    {
      // extract the value of the key and check its type
      PyObject *value = PyDict_GetItem((PyObject *)settings, key);
      if (PyList_Check(value))
      {
        // it is a list
        int listNEntries = PyList_Size(value);

        // do nothing if it is an empty list
        if (listNEntries == 0)
          return;

        // get the first value from the list
        int value = PythonUtility::getOptionListBegin<int>(settings, keyString);

        // loop over other values
        for (;
            !PythonUtility::getOptionListEnd(settings, keyString);
            PythonUtility::getOptionListNext<int>(settings, keyString, value))
        {
          values.push_back(value);
        }
      }
      else
      {
        // not a list, but a different entry (only 1 entry)
        int value = PythonUtility::getOptionInt(settings, keyString, 0);
        values.push_back(value);
      }
    }
    else
    {
      // this is no warning
      LOG(DEBUG) << "Key \"" <<keyString<< "\" not found in dict in config file. Assuming " << values;
    }
    Py_CLEAR(key);
  }
}

void PythonUtility::getOptionVector(const PyObject *settings, std::string keyString, std::vector<std::string> &values)
{
  if (settings)
  {
    // start critical section for python API calls
    PythonUtility::GlobalInterpreterLock lock;

    // check if input dictionary contains the key
    PyObject *key = PyUnicode_FromString(keyString.c_str());
    if(PyDict_Contains((PyObject *)settings, key))
    {
      // extract the value of the key and check its type
      PyObject *value = PyDict_GetItem((PyObject *)settings, key);
      if (PyList_Check(value))
      {
        // it is a list
        int listNEntries = PyList_Size(value);

        // do nothing if it is an empty list
        if (listNEntries == 0)
          return;

        // get the first value from the list
        std::string value = PythonUtility::getOptionListBegin<std::string>(settings, keyString);

        // loop over other values
        for (;
            !PythonUtility::getOptionListEnd(settings, keyString);
            PythonUtility::getOptionListNext<std::string>(settings, keyString, value))
        {
          values.push_back(value);
        }
      }
      else
      {
        // not a list, but a different entry (only 1 entry)
        std::string value = PythonUtility::getOptionString(settings, keyString, 0);
        values.push_back(value);
      }
    }
    else
    {
      // this is no warning
      LOG(DEBUG) << "Key \"" <<keyString<< "\" not found in dict in config file. Assuming " << values;
    }
    Py_CLEAR(key);
  }
}

void PythonUtility::getOptionVector(const PyObject *settings, std::string keyString, std::vector<double> &values)
{
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
        int listNEntries = PyList_Size(value);

        // do nothing if it is an empty list
        if (listNEntries == 0)
          return;

        // get the first value from the list
        double value = PythonUtility::getOptionListBegin<double>(settings, keyString);

        // loop over other values
        for (;
            !PythonUtility::getOptionListEnd(settings, keyString);
            PythonUtility::getOptionListNext<double>(settings, keyString, value))
        {
          values.push_back(value);
        }
      }
      else
      {
        // not a list, but a different entry (only 1 entry)
        double value = PythonUtility::getOptionDouble(settings, keyString, 0.0);
        values.push_back(value);
      }
    }
    else
    {
      // this is no warning
      LOG(DEBUG) << "Key \"" <<keyString<< "\" not found in dict in config file. Assuming " << values;
    }
    Py_CLEAR(key);
  }
}

void PythonUtility::getOptionVector(const PyObject *settings, std::string keyString, std::vector<PyObject *> &values)
{
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
        int listNEntries = PyList_Size(value);

        // do nothing if it is an empty list
        if (listNEntries == 0)
          return;

        // get the first value from the list
        PyObject *item = PythonUtility::getOptionListBegin<PyObject *>(settings, keyString);

        // loop over other values
        for (;
            !PythonUtility::getOptionListEnd(settings, keyString);
            PythonUtility::getOptionListNext<PyObject *>(settings, keyString, item))
        {
          values.push_back(item);
        }
      }
      else
      {
        // not a list, but a single entry (only 1 entry)
        PyObject *item = (PyObject *)settings;
        values.push_back(item);
      }
    }
    else
    {
      // this is no warning
      LOG(DEBUG) << "Key \"" <<keyString<< "\" not found in dict in config file. Assuming vector " << values;
    }
    Py_CLEAR(key);
  }
}

void PythonUtility::checkForError()
{
  if (PyErr_Occurred())
  {
    PyErr_Print();
    PyErr_Clear();
  }
}

PyObject *PythonUtility::convertToPythonList(std::vector<double> &data)
{
  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
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
  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;

  PyObject *result = PyList_New((Py_ssize_t)data.size());
  for (unsigned int i=0; i<data.size(); i++)
  {
    PyObject *item = PyLong_FromLong(data[i]);
    PyList_SetItem(result, (Py_ssize_t)i, item);    // steals reference to item
  }
  return result;    // return value: new reference
}

PyObject *PythonUtility::convertToPythonList(std::vector<global_no_t> &data)
{
  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;

  PyObject *result = PyList_New((Py_ssize_t)data.size());
  for (unsigned int i=0; i<data.size(); i++)
  {
    PyObject *item = PyLong_FromLong((long)data[i]);
    PyList_SetItem(result, (Py_ssize_t)i, item);    // steals reference to item
  }
  return result;    // return value: new reference
}

PyObject *PythonUtility::convertToPythonList(std::vector<int> &data)
{
  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;

  PyObject *result = PyList_New((Py_ssize_t)data.size());
  for (unsigned int i=0; i<data.size(); i++)
  {
    PyObject *item = PyLong_FromLong(data[i]);
    PyList_SetItem(result, (Py_ssize_t)i, item);    // steals reference to item
  }
  return result;    // return value: new reference
}

PyObject *PythonUtility::convertToPythonList(std::vector<bool> &data)
{
  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;

  PyObject *result = PyList_New((Py_ssize_t)data.size());
  for (unsigned int i=0; i<data.size(); i++)
  {
    PyObject *item = (data[i]? Py_True : Py_False);
    PyList_SetItem(result, (Py_ssize_t)i, item);    // steals reference to item
  }
  return result;    // return value: new reference
}

PyObject *PythonUtility::convertToPythonList(unsigned int nEntries, double* data)
{
  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
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
  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
#if PY_MAJOR_VERSION >= 3
  PyObject *asciiString = PyUnicode_AsASCIIString(object);
  std::string result = PyBytes_AsString(asciiString);
  Py_DECREF(asciiString);
#else
  std::string result = pyUnicodeToString(object);
#endif

  return result;
}


// python GIL handling is only needed for multi-threading. It does not work properly with openmp.
//int PythonUtility::GlobalInterpreterLock::nGILS_ = 0;
/*std::recursive_mutex PythonUtility::GlobalInterpreterLock::mutex_;
std::unique_lock<std::recursive_mutex> PythonUtility::GlobalInterpreterLock::lock_(
  PythonUtility::GlobalInterpreterLock::mutex_, std::defer_lock);
*/
/*
omp_nest_lock_t PythonUtility::GlobalInterpreterLock::lock_;
bool PythonUtility::GlobalInterpreterLock::lockInitialized_ = false;
std::map<int, int> PythonUtility::GlobalInterpreterLock::nGilsThreads_;
*/
PythonUtility::GlobalInterpreterLock::GlobalInterpreterLock()
{
  // start critical section for python interpreter
  // store GlobalInterpreterLock state
  //nGILS_++;
  //LOG(INFO) << omp_get_thread_num() << ": nGILS: " << nGILS_;
  //lock_.lock();
  /*if (nGilsThreads_.find(omp_get_thread_num()) == nGilsThreads_.end())
  {
    nGilsThreads_[omp_get_thread_num()] = 0;
  }
  
  nGilsThreads_[omp_get_thread_num()]++;
  
  if (!lockInitialized_)
  {
    omp_init_nest_lock(&lock_);
    lockInitialized_ = true;
  }
  
  omp_set_nest_lock(&lock_);
  */
  //if (nGilsThreads_[omp_get_thread_num()] == 1)
  
#if 0  
  LOG(INFO) << "[" << omp_get_thread_num() << "/" << omp_get_thread_num() << "] (wait for GIL)";
  gstate_ = PyGILState_Ensure();
  
  LOG(INFO) << "[" << omp_get_thread_num() << "/" << omp_get_thread_num() << "] (acquired)";
#endif  
  
  //mainThreadState_ = PyEval_SaveThread();

  //Py_BEGIN_ALLOW_THREADS
  //LOG(DEBUG) << "PyGILState_Ensure";
}

PythonUtility::GlobalInterpreterLock::~GlobalInterpreterLock()
{
 /* 
  nGilsThreads_[omp_get_thread_num()]--;
  
  // Release the thread. No Python API allowed beyond this point.
  if (nGilsThreads_[omp_get_thread_num()] == 0)
  {*/
#if 0 
  PyGILState_Release(gstate_);
#endif
  /*}
  nGILS_--;
  //lock_.unlock();
  
  omp_unset_nest_lock(&lock_);
  
  //LOG(INFO) << omp_get_thread_num() << ": nGILS: " << nGILS_ << "(released)";
  
  //PyEval_RestoreThread(mainThreadState_); 
  //Py_END_ALLOW_THREADS
  //LOG(DEBUG) << "PyGILState_Release";
  */
}

std::ostream &operator<<(std::ostream &stream, PyObject *object)
{
  stream << PythonUtility::getString(object);
  return stream;
}
