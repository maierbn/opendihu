#include "control/python_utility.h"

#include <fstream>
#include <iostream>
#include <iomanip>

#include "Python.h"
#include "easylogging++.h"

PyObject *PythonUtility::itemList = NULL;
int PythonUtility::itemListIndex = 0;
PyObject *PythonUtility::list = NULL;
int PythonUtility::listIndex = 0;

template<>
int PythonUtility::convertFromPython(PyObject *object, int defaultValue)
{
  if (PyInt_Check(object))
  {
    long valueLong = PyInt_AsLong(object);
    return int(valueLong);
  }
  else if (PyFloat_Check(object))
  {
    double valueDouble = PyFloat_AsDouble(object);
      
    if (double(int(valueDouble)) != valueDouble)      // if value is not e.g. 2.0 
    {
      LOG(WARNING) << "convertFromPython: object is no int.";
    }
    
    return int(valueDouble);
  }
  else if (PyString_Check(object))
  {
    std::string valueString = PyString_AsString(object);
    return atoi(valueString.c_str());
  }
  else
  {
    LOG(WARNING) << "convertFromPython: object is no int.";
  }
  return defaultValue;
}

template<>
double PythonUtility::convertFromPython(PyObject *object, double defaultValue)
{
  if (PyFloat_Check(object))
  {
    double valueDouble = PyFloat_AsDouble(object);
      
    return valueDouble;
  }
  else if (PyInt_Check(object))
  {
    long valueLong = PyInt_AsLong(object);
    return double(valueLong);
  }
  else if (PyString_Check(object))
  {
    std::string valueString = PyString_AsString(object);
    return atof(valueString.c_str());
  }
  else
  {
    LOG(WARNING) << "convertFromPython: object is no double.";
  }
  return defaultValue;
}

template<>
std::string PythonUtility::convertFromPython(PyObject *object, std::string defaultValue)
{
  if (PyString_Check(object))
  {
    std::string valueString = PyString_AsString(object);
    return valueString;
  }
  else
  {
    LOG(WARNING) << "convertFromPython: object is no std::string.";
  }
  return defaultValue;
}

template<>
int PythonUtility::convertFromPython(PyObject *object)
{
  return convertFromPython<int>(object, 0);
}

template<>
bool PythonUtility::convertFromPython(PyObject *object, bool defaultValue)
{
  if (PyBool_Check(object))
  {
    if (object == Py_True)
    {
      return true;
    }
    return false;
  }
  if (PyInt_Check(object))
  {
    long valueLong = PyInt_AsLong(object);
    return bool(valueLong);
  }
  else if (PyString_Check(object))
  {
    std::string valueString = PyString_AsString(object); 
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
    LOG(WARNING) << "convertFromPython: object is no bool.";
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

PyObject *PythonUtility::extractDict(PyObject *dict, std::string keyString)
{
  if (dict)
  {
    // check if input dictionary contains the key
    PyObject *key = PyString_FromString(keyString.c_str());
    
    if(PyDict_Contains(dict, key))
    {
      return PyDict_GetItem(dict, key);
    }
    else
    {
      LOG(WARNING) << "Dict does not contain Key \""<<keyString<<"\"!"<<std::endl;
      return NULL;
    }
  }
  return NULL;
}

double PythonUtility::getOptionDouble(PyObject* settings, std::string keyString, double defaultValue, ValidityCriterion validityCriterion)
{
  double result = defaultValue;
  
  if (!settings)
    return result;
  
  // check if input dictionary contains the key
  PyObject *key = PyString_FromString(keyString.c_str());
  if(PyDict_Contains(settings, key))
  {
    // extract the value of the key and check its type
    PyObject *value = PyDict_GetItem(settings, key);
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
  case None:
    break;
  };
  
  return result;
}

int PythonUtility::getOptionInt(PyObject *settings, std::string keyString, int defaultValue, ValidityCriterion validityCriterion)
{
  int result = defaultValue;
  
  if (!settings)
    return result;
  
  // check if input dictionary contains the key
  PyObject *key = PyString_FromString(keyString.c_str());
  if(PyDict_Contains(settings, key))
  { 
    // extract the value of the key and check its type
    PyObject *value = PyDict_GetItem(settings, key);
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
      if (result <= 0.0)
      {
        LOG(WARNING)<<"value "<<result<<" of Key \""<<keyString<<"\" is invalid (not positive). Using default value "
          <<defaultValue<<".";
        result = defaultValue;
      }
          break;
    case None:
      break;
  };
  
  return result;
}

bool PythonUtility::getOptionBool(PyObject *settings, std::string keyString, bool defaultValue)
{
  int result = defaultValue;
  
  if (!settings)
    return result;
  
  // check if input dictionary contains the key
  PyObject *key = PyString_FromString(keyString.c_str());
  if(PyDict_Contains(settings, key))
  { 
    // extract the value of the key and check its type
    PyObject *value = PyDict_GetItem(settings, key);
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
  return result;
}

std::string PythonUtility::getOptionString(PyObject *settings, std::string keyString, std::string defaultValue)
{
  std::string result = defaultValue;
  
  if (!settings)
    return result;
  
  // check if input dictionary contains the key
  PyObject *key = PyString_FromString(keyString.c_str());
  if(PyDict_Contains(settings, key))
  { 
    // extract the value of the key and check its type
    PyObject *value = PyDict_GetItem(settings, key);
    
    // convert to std::string or try default value
    result = convertFromPython<std::string>(value, defaultValue);
  }
  else
  {
    LOG(WARNING)<<"Key \""<<keyString<<"\" not found in dict in config file.";
  }
  
  return result;
}

void PythonUtility::printDict(PyObject *dict, int indent)
{
  // iterate over top level key-value pairs
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  
  while (PyDict_Next(dict, &pos, &key, &value))
  {
    LOG(INFO)<<std::string(indent, ' ');
    
    if (!PyString_Check(key))
    {
      LOG(INFO)<<"key is not a string"<<std::endl;
    }
    else
    {
      std::string keyString = PyString_AsString(key);
      LOG(INFO)<<keyString<<": ";
    }
                
    if (PyString_CheckExact(value))
    {
      std::string valueString = PyString_AsString(value);
      LOG(INFO)<<"\""<<valueString<<"\""<<std::endl;
    }
    else if (PyInt_CheckExact(value))
    {
      long valueLong = PyInt_AsLong(value);
      LOG(INFO)<<valueLong<<std::endl;
    }
    else if (PyLong_CheckExact(value))
    {
      long valueLong = PyLong_AsLong(value);
      LOG(INFO)<<valueLong<<std::endl;
    }
    else if (PyFloat_CheckExact(value))
    {
      double valueDouble = PyFloat_AsDouble(value);
      LOG(INFO)<<valueDouble<<std::endl;
    }
    else if (PyBool_Check(value))
    {
      bool valueBool = PyObject_IsTrue(value);
      LOG(INFO)<<std::boolalpha<<valueBool<<std::endl;
    }
    else if(PyDict_CheckExact(value))
    {
      LOG(INFO)<<"{"<<std::endl;
      printDict(value, indent+2);
      LOG(INFO)<<std::string(indent, ' ')<<"}"<<std::endl;
    }
    else
    {
      LOG(INFO)<<"<unknown>"<<std::endl;
    }
  }
}

bool PythonUtility::getOptionDictEnd(PyObject *settings, std::string keyString)
{
  if (!itemList)
    return true;
  return itemListIndex >= PyList_Size(itemList);
}

bool PythonUtility::getOptionListEnd(PyObject *settings, std::string keyString)
{
  if (!list)
    return true;
  return listIndex >= PyList_Size(list);
}

void PythonUtility::getOptionVector(PyObject* settings, std::string keyString, int nEntries, std::vector<double> &values)
{
  values.resize(nEntries);
  
  // check if input dictionary contains the key
  PyObject *key = PyString_FromString(keyString.c_str());
  if(PyDict_Contains(settings, key))
  { 
    // extract the value of the key and check its type
    PyObject *value = PyDict_GetItem(settings, key);
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
}