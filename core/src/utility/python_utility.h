#pragma once

#include <Python.h>  // has to be the first included header

#include <iostream>
#include <array>
#include <vector>

// With python3+ PyString_* was renamed to PyBytes_*
//(This ugly check should be removed when decided if python2.7 or python3 will be used. Do that now and remove it now!)
#if PY_MAJOR_VERSION >= 3
#define PyString_Check PyBytes_Check
#define PyString_CheckExact PyBytes_CheckExact
#define PyString_FromString PyBytes_FromString
#define PyString_AsString PyBytes_AsString
#endif

/** Utility class that handles parsing of python config data to c type objects
 */
class PythonUtility
{
public:
  enum ValidityCriterion
  {
    None,
    Positive,
    NonNegative,
    Between1And3
  };
 
  //! constructor, initialize context, parse command line parameters and input file
  PythonUtility();
    
  //! checks if the settings contain the given key, no warning is printed
  static bool containsKey(const PyObject *settings, std::string key);
  
  //! given a python dictionary in settings, extract the value of given key and check if it is again a dict. Returns NULL, if the key does not exist. Then also a warning is printed.
  static PyObject *getOptionPyObject(const PyObject *settings, std::string key);
  
  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  static double getOptionDouble(const PyObject *settings, std::string key, double defaultValue, ValidityCriterion validityCriterion = None);
  
  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  static int getOptionInt(const PyObject *settings, std::string key, int defaultValue, ValidityCriterion validityCriterion = None);
  
  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  static bool getOptionBool(const PyObject *settings, std::string key, bool defaultValue=true);
  
  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  static std::string getOptionString(const PyObject *settings, std::string key, std::string defaultValue);
  
  //! return the option value given by key in the python dictionary settings. If not found, return NULL
  static PyObject *getOptionFunction(const PyObject *settings, std::string key);
  
  //! return the option value as array given by key in the python dictionary settings. If not found, return the defaultValue, also check if validityCriterion is met
  template<class ValueType, int D>
  static std::array<ValueType, D> getOptionArray(PyObject* settings, std::string keyString, std::array<ValueType, D> defaultValue, 
                                                ValidityCriterion validityCriterion = None);
  
  //! return the option value as array given by key in the python dictionary settings. If not found, return the defaultValue, also check if validityCriterion is met
  template<class ValueType, int D>
  static std::array<ValueType, D> getOptionArray(PyObject* settings, std::string keyString, ValueType defaultValue, 
                                                ValidityCriterion validityCriterion = None);
  
  //! Consider a Python dictionary in settings with key keyString. Set an internal iterator to the beginning and return the first key,value pair
  template<typename Key, typename Value>
  static std::pair<Key, Value> getOptionDictBegin(const PyObject *settings, std::string keyString);
  
  //! If the internal iterator on the current dictionary is at the end of the dictionary
  static bool getOptionDictEnd(const PyObject *settings, std::string keyString);
  
  //! Increment the internal iterator of the dictionary and set the next key,value pair in nextPair
  template<typename Key, typename Value>
  static void getOptionDictNext(const PyObject *settings, std::string keyString, std::pair<Key, Value> &nextPair);
  
  //! Consider a Python list in settings with key keyString. Set an internal iterator to the beginning and return the first key,value pair
  template<typename Value>
  static Value getOptionListBegin(const PyObject *settings, std::string keyString);
  
  //! If the internal iterator on the current list is at the end of the list
  static bool getOptionListEnd(const PyObject *settings, std::string keyString);
  
  //! Increment the internal iterator of the list and set the next key,value pair in nextPair
  template<typename Value>
  static void getOptionListNext(const PyObject *settings, std::string keyString, Value &value);

  //! extract a vector with exactly the specified number of nEntries, can be a dict or list, not specified entries are set to 0
  static void getOptionVector(const PyObject *settings, std::string keyString, int nEntries, std::vector<double> &values);
  
  //! print python dictionary to stdout
  static void printDict(PyObject *dict, int indent=0);
  
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
  template<typename T>
  static T convertFromPython(PyObject *object, T defaultValue);
  
  //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
  template<typename T>
  static T convertFromPython(PyObject *object);
  
  //! create a python list out of the double vector
  static PyObject *convertToPythonList(std::vector<double> &data);
  
  //! create a python list out of the long vector
  static PyObject *convertToPythonList(std::vector<long> &data);
  
  //! create a python list out of the long array
  template<int D>
  static PyObject *convertToPythonList(std::array<long,D> &data);
  
  //! create a python list out of the long vector
  static PyObject *convertToPythonList(unsigned int nEntries, double *data);
 
  //! convert a python list to a std::array, use default value when python list is shorter than the resulting array
  template<class ValueType, int D>
  static std::array<ValueType, D> convertFromPython(PyObject *object, ValueType defaultValue);
  
  //! convert a python list to a std::array with specified default values if the python list is shorter than the array
  template<class ValueType, int D>
  static std::array<ValueType, D> convertFromPython(PyObject *object, std::array<ValueType, D> defaultValue);
  
  //! convert a python list to a std::array
  template<class ValueType, int D>
  static std::array<ValueType, D> convertFromPython(PyObject *object);
  
private:
  
  static PyObject *itemList;    ///< list of items (key,value) for dictionary,  to use for getOptionDictBegin, getOptionDictEnd, getOptionDictNext
  static int itemListIndex;     ///< current index of itemList
  
  static PyObject *list;      ///< python list to use for getOptionListBegin, getOptionListEnd, getOptionListNext
  static int listIndex;       ///< current index for list
};

#include "utility/python_utility.tpp"
