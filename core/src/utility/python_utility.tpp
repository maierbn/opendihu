#include "utility/python_utility.h"

#include <Python.h>
#include "easylogging++.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>

#include "utility/vector_operators.h"
#include "control/python_config/settings_file_name.h"

template<typename Key, typename Value>
std::pair<Key, Value> PythonUtility::getOptionDictBegin(const PyObject *settings, std::string keyString, std::string pathString)
{
  std::pair<Key, Value> firstEntry;

  if (settings && PyDict_Check(settings))
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;
  
    // check if input dictionary contains the key
    PyObject *key = PyUnicode_FromString(keyString.c_str());
    if (PyDict_Contains((PyObject *)settings, key))
    {
      //PythonUtility::printDict((PyObject *)settings);

      PyObject *dict = PyDict_GetItem((PyObject *)settings, key);
      if (PyDict_Check(dict))
      {
        Py_CLEAR(key);
        Py_CLEAR(itemList);
        itemList = PyDict_Items(dict);
        itemListIndex = 0;

        if (PyList_Check(itemList))
        {

          if (itemListIndex < PyList_Size(itemList))
          {
            PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
            PyObject *tuple_key = PyTuple_GetItem(tuple, (Py_ssize_t)0);
            PyObject *tuple_value = PyTuple_GetItem(tuple, (Py_ssize_t)1);

            firstEntry = std::pair<Key, Value>(convertFromPython<Key>::get(tuple_key), convertFromPython<Value>::get(tuple_value));
            return firstEntry;
          }
        }
        else
        {
          LOG(WARNING) << pathString << "[\"" << keyString << "\"] is not a dict";
        }
      }
      else
      {
        LOG(WARNING) << "Entry " << pathString << "[\"" << keyString << "\"] is not a dict.";
      }
    }
    else
    {
      LOG(WARNING) << pathString << "[\"" << keyString << "\"] not set in \"" << Control::settingsFileName << "\"" << std::endl;
    }
  }

  return firstEntry;
}

template<typename Key, typename Value>
void PythonUtility::getOptionDictNext(const PyObject *settings, std::string keyString, std::string pathString, std::pair<Key, Value> &nextPair)
{
  itemListIndex++;

  // start critical section for python API calls
  // PythonUtility::GlobalInterpreterLock lock;
  
  if (itemListIndex < PyList_Size(itemList))
  {
    PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
    PyObject *key = PyTuple_GetItem(tuple, (Py_ssize_t)0);
    PyObject *value = PyTuple_GetItem(tuple, (Py_ssize_t)1);

    nextPair = std::pair<Key, Value>(convertFromPython<Key>::get(key), convertFromPython<Value>::get(value));

  }
}

template<typename Value>
Value PythonUtility::getOptionListBegin(const PyObject *settings, std::string keyString, std::string pathString)
{
  if (settings && PyDict_Check(settings))
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;
  
    // check if input dictionary contains the key
    PyObject *key = PyUnicode_FromString(keyString.c_str());
    if (PyDict_Contains((PyObject *)settings, key))
    {
      // check if it is a list
      list = PyDict_GetItem((PyObject *)settings, key);
      if (PyList_Check(list))
      {
        listIndex = 0;

        if (listIndex < PyList_Size(list))
        {
          PyObject *item = PyList_GetItem(list, (Py_ssize_t)listIndex);

          Py_CLEAR(key);
          return convertFromPython<Value>::get(item);
        }
      }
      else
      {
        LOG(WARNING) << "" << pathString << "[\"" << keyString << "\"] is not a list!";
        Py_CLEAR(key);
        return convertFromPython<Value>::get(list);
      }
    }
    else
    {
      LOG(WARNING) << pathString << "[\"" << keyString << "\"] not found in config file.";
    }

    Py_CLEAR(key);
  }

  return Value();
}

template<typename Value>
void PythonUtility::getOptionListNext(const PyObject *settings, std::string keyString, std::string pathString, Value &value)
{
  listIndex++;

  // start critical section for python API calls
  // PythonUtility::GlobalInterpreterLock lock;
  
  if (listIndex < PyList_Size(list))
  {
    PyObject *item = PyList_GetItem(list, (Py_ssize_t)listIndex);

    value = convertFromPython<Value>::get(item);
  }
}

template<class ValueType, int D>
std::array<ValueType, D> PythonUtility::getOptionArray(PyObject* settings, std::string keyString, std::string pathString,
                                                      ValueType defaultValue, ValidityCriterion validityCriterion)
{
  std::array<ValueType,(int)D> defaultValueArray = {};
  defaultValueArray.fill(defaultValue);
  return PythonUtility::getOptionArray<ValueType,D>(settings, keyString, pathString, defaultValueArray, validityCriterion);
}

template<class ValueType, int D>
std::array<ValueType, D> PythonUtility::getOptionArray(PyObject* settings, std::string keyString, std::string pathString,
                                                      std::array<ValueType, D> defaultValue, ValidityCriterion validityCriterion)
{
  std::array<ValueType, D> result = defaultValue;

  if (!settings || !PyDict_Check(settings))
  {
    return result;
  }

  if (settings)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;
  
    // check if input dictionary contains the key
    PyObject *key = PyUnicode_FromString(keyString.c_str());
    if (PyDict_Contains((PyObject *)settings, key))
    {
      // extract the value of the key and check its type
      PyObject *value = PyDict_GetItem((PyObject *)settings, key);
      result = PythonUtility::convertFromPython<std::array<ValueType,D>>::get(value, defaultValue);
    }
    else
    {
      LOG(WARNING) << pathString << "[\"" << keyString << "\"] not found in config, assuming default values " << defaultValue << ".";

      Py_CLEAR(key);
      return defaultValue;
    }

    Py_CLEAR(key);
  }

  switch(validityCriterion)
  {
    case PythonUtility::Positive:
      for (int i=0; i<D; i++)
      {
       if (result[i] <= 0.0)
       {
         LOG(WARNING) << "Value " <<result[i]<< " of " << pathString << "[\"" << keyString << "\"] is invalid (not positive). Using default value "
           << defaultValue[i]<< ".";
         result[i] = defaultValue[i];
       }
      }
    case PythonUtility::NonNegative:
      for (int i=0; i<D; i++)
      {
       if (result[i] < 0.0)
       {
         LOG(WARNING) << "Value " <<result[i]<< " of " << pathString << "[\"" << keyString << "\"] is invalid (not non-negative). Using default value "
           << defaultValue[i]<< ".";
         result[i] = defaultValue[i];
       }
      }

    break;
    case PythonUtility::Between1And3:
      for (int i=0; i<D; i++)
      {
       if (result[i] < 1.0 || result[i] > 3.0)
       {
         LOG(WARNING) << "Value " <<result[i]<< " of " << pathString << "[\"" << keyString << "\"] is invalid (not between 1 and 3). Using default value "
           << defaultValue[i]<< ".";
         result[i] = defaultValue[i];
       }
      }

    break;
    case PythonUtility::None:
      break;
  };

  return result;
}

template<int D>
PyObject *PythonUtility::convertToPythonList(std::array<long,D> &data)
{
  // start critical section for python API calls
  // PythonUtility::GlobalInterpreterLock lock;

  PyObject *result = PyList_New((Py_ssize_t)D);
  for (unsigned int i=0; i<D; i++)
  {
    PyObject *item = PyLong_FromLong(data[i]);
    PyList_SetItem(result, (Py_ssize_t)i, item);    // steals reference to item
  }
  return result;    // return value: new reference
}

template<int D>
PyObject *PythonUtility::convertToPythonList(std::array<bool,D> &data)
{
  // start critical section for python API calls
  // PythonUtility::GlobalInterpreterLock lock;

  PyObject *result = PyList_New((Py_ssize_t)D);
  for (unsigned int i=0; i<D; i++)
  {
    PyObject *item = (data[i]? Py_True : Py_False);
    PyList_SetItem(result, (Py_ssize_t)i, item);    // steals reference to item
  }
  return result;    // return value: new reference
}
