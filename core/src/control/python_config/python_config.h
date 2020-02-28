#pragma once

#include <Python.h>  // has to be the first included header
#include "utility/python_utility.h"

class PythonConfig
{
public:
  //! constructor from python object
  PythonConfig(PyObject *specificSettings = NULL);

  //! constructor as sub scope of another python config
  PythonConfig(const PythonConfig &rhs, std::string key);

  //! constructor as sub scope of another python config which is a list
  PythonConfig(const PythonConfig &rhs, int i);

  //! constructor directly from PyObject*, path from rhs + key
  PythonConfig(const PythonConfig &rhs, std::string key, PyObject *config);

  //! constructor directly from PyObject*, path from rhs + key + key2
  PythonConfig(const PythonConfig &rhs, std::string key1, std::string key2, PyObject *config);

  //! default copy constructor
  PythonConfig(PythonConfig const &) = default;

  //! default move constructor
  PythonConfig(PythonConfig&&) = default;

  //! destructor
  ~PythonConfig();

  //! default assignment operator
  PythonConfig &operator=(const PythonConfig &);

  //! return the actual PyObject object
  PyObject *pyObject() const;

  //! set the internal pyObject object
  void setPyObject(PyObject *pyObject);

  //! iterator to beginning of path vector, which contains the keys from the top-level config object to the current context
  std::vector<std::string>::const_iterator pathBegin() const;

  //! iterator to end of path vector, which contains the keys from the top-level config object to the current context
  std::vector<std::string>::const_iterator pathEnd() const;

  //! checks if this settings contain the given key, no warning is printed (hasOption)
  bool hasKey(std::string key) const;

  //! return all keys of the current dict as vector of strings
  void getKeys(std::vector<std::string> &keys);

  //! given a python dictionary in settings, extract the value of given key and check if it is again a dict. Returns NULL, if the key does not exist. Then also a warning is printed.
  PyObject *getOptionPyObject(std::string key, PyObject *defaultValue = Py_None) const;

  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  double getOptionDouble(std::string key, double defaultValue, PythonUtility::ValidityCriterion validityCriterion = PythonUtility::ValidityCriterion::None) const;

  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  int getOptionInt(std::string key, int defaultValue, PythonUtility::ValidityCriterion validityCriterion = PythonUtility::ValidityCriterion::None) const;

  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  bool getOptionBool(std::string key, bool defaultValue = true) const;

  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  std::string getOptionString(std::string key, std::string defaultValue) const;

  //! return the option value given by key in the python dictionary settings. If not found, return NULL
  PyObject *getOptionFunction(std::string key) const;

  //! return the option value as array given by key in the python dictionary settings. If not found, return the defaultValue, also check if validityCriterion is met
  template<class ValueType, int D>
  std::array<ValueType, D> getOptionArray(std::string keyString, std::array<ValueType, D> defaultValue,
                                          PythonUtility::ValidityCriterion validityCriterion = PythonUtility::ValidityCriterion::None) const;

  //! return the option value as array given by key in the python dictionary settings. If not found, return the defaultValue, also check if validityCriterion is met
  template<class ValueType, int D>
  std::array<ValueType, D> getOptionArray(std::string keyString, ValueType defaultValue,
                                          PythonUtility::ValidityCriterion validityCriterion = PythonUtility::ValidityCriterion::None) const;

  //! Consider a Python dictionary in settings with key keyString. Set an internal iterator to the beginning and return the first key,value pair
  template<typename Key, typename Value>
  std::pair<Key, Value> getOptionDictBegin(std::string keyString) const;

  //! If the internal iterator on the current dictionary is at the end of the dictionary
  bool getOptionDictEnd(std::string keyString) const;

  //! Increment the internal iterator of the dictionary and set the next key,value pair in nextPair
  template<typename Key, typename Value>
  void getOptionDictNext(std::string keyString, std::pair<Key, Value> &nextPair) const;

  //! Consider a Python list in settings with key keyString. Set an internal iterator to the beginning and return the first key,value pair
  template<typename Value>
  Value getOptionListBegin(std::string keyString) const;

  //! If the internal iterator on the current list is at the end of the list
  bool getOptionListEnd(std::string keyString) const;

  //! Increment the internal iterator of the list and set the next key,value pair in nextPair
  template<typename Value>
  void getOptionListNext(std::string keyString, Value &value) const;

  //! extract a vector with exactly the specified number of nEntries, can be a dict or list, not specified entries are set to 0
  void getOptionVector(std::string keyString, int nEntries, std::vector<double> &values) const;

  //! extract a vector with exactly the specified number of nEntries, must be a list, not specified entries are set to {0}
  template<typename T>
  void getOptionVector(std::string keyString, int nEntries, std::vector<T> &values) const;

  //! extract a vector with unknown number of nEntries, must be a list
  template<typename T>
  void getOptionVector(std::string keyString, std::vector<T> &values) const;

  //! extract a vector with unknown number of nEntries, must be a list
  void getOptionVector(std::string keyString, std::vector<PythonConfig> &values) const;

  //! get the python code of this object, e.g. config["key1"]["key2"] etc.
  std::string getStringPath() const;

protected:

  PyObject *pythonConfig_;    ///< the python config dictionary of the current context (i.e. may be a sub-dict of the global config)
  std::vector<std::string> path_;   ///< the key words of the python config down to the current scope
};

//! output operator
std::ostream &operator<<(std::ostream &stream, const PythonConfig rhs);

#include "control/python_config/python_config.tpp"
