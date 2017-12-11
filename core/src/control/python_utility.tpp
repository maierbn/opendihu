#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>

#include "Python.h"
#include "easylogging++.h"

template<typename Key, typename Value>
std::pair<Key, Value> PythonUtility::getOptionDictBegin(const PyObject *settings, std::string keyString)
{
  //PyDict_Check(settings
  
  std::pair<Key, Value> firstEntry;
 
  if (settings)
  {
    // check if input dictionary contains the key
    PyObject *key = PyString_FromString(keyString.c_str());
    if(PyDict_Contains((PyObject *)settings, key))
    {
      PyObject *dict = PyDict_GetItem((PyObject *)settings, key);
      itemList = PyDict_Items(dict);
      itemListIndex = 0;
      
      if (itemListIndex < PyList_Size(itemList))
      {
        PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
        PyObject *key = PyTuple_GetItem(tuple, (Py_ssize_t)0);
        PyObject *value = PyTuple_GetItem(tuple, (Py_ssize_t)1);
        
        return std::pair<Key, Value>(convertFromPython<Key>(key), convertFromPython<Value>(value));
      }
    }
    else
    {
      LOG(WARNING)<<"Warning: key \""<<keyString<<"\" not found in dict in config file"<<std::endl;
    }
  }
  
  return firstEntry;
}

template<typename Key, typename Value>
void PythonUtility::getOptionDictNext(const PyObject *settings, std::string keyString, std::pair<Key, Value> &nextPair)
{
  itemListIndex++;
  
  if (itemListIndex < PyList_Size(itemList))
  {
    PyObject *tuple = PyList_GetItem(itemList, (Py_ssize_t)itemListIndex);
    PyObject *key = PyTuple_GetItem(tuple, (Py_ssize_t)0);
    PyObject *value = PyTuple_GetItem(tuple, (Py_ssize_t)1);
    
    nextPair = std::pair<Key, Value>(convertFromPython<Key>(key), convertFromPython<Value>(value));
  }
}

template<typename Value>
Value PythonUtility::getOptionListBegin(const PyObject *settings, std::string keyString)
{
  if (settings)
  {
    // check if input dictionary contains the key
    PyObject *key = PyString_FromString(keyString.c_str());
    if(PyDict_Contains((PyObject *)settings, key))
    {
      // check if it is a list
      list = PyDict_GetItem((PyObject *)settings, key);
      if (PyList_Check(list))
      {
        listIndex = 0;
        
        if (listIndex < PyList_Size(list))
        {
          PyObject *item = PyList_GetItem(list, (Py_ssize_t)listIndex);
          
          return convertFromPython<Value>(item);
        }
      }
      else
      {
        LOG(WARNING)<<"Key \""<<keyString<<"\" is not a list!";
        return convertFromPython<Value>(list);
      }
    }
    else
    {
      LOG(WARNING)<<"Warning: key \""<<keyString<<"\" not found in dict in config file"<<std::endl;
    }
  }
  
  return Value();
}

template<typename Value>
void PythonUtility::getOptionListNext(const PyObject *settings, std::string keyString, Value &value)
{
  listIndex++;
  
  if (listIndex < PyList_Size(list))
  {
    PyObject *item = PyList_GetItem(list, (Py_ssize_t)listIndex);
    
    value = convertFromPython<Value>(item);
  }
}

template<class ValueType, unsigned long D>
std::array<ValueType, D> PythonUtility::convertFromPython(PyObject *object, std::array<ValueType, D> defaultValue)
{
  std::array<ValueType, D> result;
  if (PyList_Check(object))
  {
    unsigned long i = 0;
    unsigned long iEnd = std::min((unsigned long)PyList_Size(object), D);
    
    for(;i < iEnd; i++)
    {
      result[i] = PythonUtility::convertFromPython<ValueType>(PyList_GetItem(object, (Py_ssize_t)i), defaultValue[i]);
    }
    
    // fill rest of values with default values
    for(;i < D; i++)
    {
      result[i] = defaultValue[i];
    }
    return result;
  }
  else
  {
    ValueType valueDouble = PythonUtility::convertFromPython<ValueType>(object, defaultValue[0]);
    
    result[0] = valueDouble;
    std::copy(defaultValue.begin()+1, defaultValue.end(), result.begin()+1);
    
    return result;
  }
  return defaultValue;
}

template<class ValueType, unsigned long D>
std::array<ValueType, D> PythonUtility::convertFromPython(PyObject *object)
{
  std::array<ValueType, D> defaultValue;
  defaultValue.fill(0.0);
  return convertFromPython<ValueType, D>(object, defaultValue);
}


template<class ValueType, unsigned long D>
std::array<ValueType, D> PythonUtility::getOptionArray(PyObject* settings, std::string keyString,
                                                      ValueType defaultValue, ValidityCriterion validityCriterion)
{ 
  std::array<ValueType, (unsigned long int)D> defaultValueArray = {};
  defaultValueArray.fill(defaultValue);
  return getOptionArray(settings, keyString, defaultValueArray, validityCriterion);
}

template<class ValueType, unsigned long D>
std::array<ValueType, D> PythonUtility::getOptionArray(PyObject* settings, std::string keyString,
                                                      std::array<ValueType, D> defaultValue, ValidityCriterion validityCriterion)
{ 
  std::array<ValueType, D> result = {};
 
  if (settings)
  { 
    // check if input dictionary contains the key
    PyObject *key = PyString_FromString(keyString.c_str());
    if(PyDict_Contains((PyObject *)settings, key))
    {
      // extract the value of the key and check its type
      PyObject *value = PyDict_GetItem((PyObject *)settings, key);
      result = PythonUtility::convertFromPython<ValueType, D>(value, defaultValue);
    }
    else
    {
      LOG(WARNING)<<"Warning: key \""<<keyString<<"\" not found in dict in config file"<<std::endl;
      return defaultValue;
    }
  }
  
  switch(validityCriterion)
  {
    case PythonUtility::Positive:
      for (unsigned long i=0; i<D; i++)
      {
       if (result[i] <= 0.0)
       {
         LOG(WARNING)<<"Warning: value "<<result[i]<<" of key \""<<keyString<<"\" is invalid (not positive). Using default value "
           <<defaultValue[i]<<".";
         result[i] = defaultValue[i];
       } 
      }
    case PythonUtility::NonNegative:
      for (unsigned long i=0; i<D; i++)
      {
       if (result[i] < 0.0)
       {
         LOG(WARNING)<<"Warning: value "<<result[i]<<" of key \""<<keyString<<"\" is invalid (not non-negative). Using default value "
           <<defaultValue[i]<<".";
         result[i] = defaultValue[i];
       } 
      }
      
    break;
    case PythonUtility::None:
      break;
  };
  
  return result;
}
