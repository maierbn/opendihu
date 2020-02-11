#include "control/python_config.h"

#include "easylogging++.h"
#include "utility/python_utility.h"

//! return the option value as array given by key in the python dictionary settings. If not found, return the defaultValue, also check if validityCriterion is met
template<class ValueType, int D>
std::array<ValueType, D> PythonConfig::
getOptionArray(std::string keyString, std::array<ValueType, D> defaultValue,
               PythonUtility::ValidityCriterion validityCriterion) const
{
  std::string pathString = this->getStringPath();
  return PythonUtility::getOptionArray<ValueType,D>(this->pythonConfig_, keyString, pathString, defaultValue, validityCriterion);
}

//! return the option value as array given by key in the python dictionary settings. If not found, return the defaultValue, also check if validityCriterion is met
template<class ValueType, int D>
std::array<ValueType, D> PythonConfig::
getOptionArray(std::string keyString, ValueType defaultValue,
              PythonUtility::ValidityCriterion validityCriterion) const
{
  std::string pathString = this->getStringPath();
  return PythonUtility::getOptionArray<ValueType,D>(this->pythonConfig_, keyString, pathString, defaultValue, validityCriterion);
}

//! Consider a Python dictionary in settings with key keyString. Set an internal iterator to the beginning and return the first key,value pair
template<typename Key, typename Value>
std::pair<Key, Value> PythonConfig::
getOptionDictBegin(std::string keyString) const
{
  std::string pathString = this->getStringPath();
  return PythonUtility::getOptionDictBegin<Key, Value>(this->pythonConfig_, keyString, pathString);
}

//! Increment the internal iterator of the dictionary and set the next key,value pair in nextPair
template<typename Key, typename Value>
void PythonConfig::
getOptionDictNext(std::string keyString, std::pair<Key, Value> &nextPair) const
{
  std::string pathString = this->getStringPath();
  PythonUtility::getOptionDictNext<Key,Value>(this->pythonConfig_, keyString, pathString, nextPair);
}

//! Consider a Python list in settings with key keyString. Set an internal iterator to the beginning and return the first key,value pair
template<typename Value>
Value PythonConfig::
getOptionListBegin(std::string keyString) const
{
  std::string pathString = this->getStringPath();
  return PythonUtility::getOptionListBegin<Value>(this->pythonConfig_, keyString, pathString);
}

//! Increment the internal iterator of the list and set the next key,value pair in nextPair
template<typename Value>
void PythonConfig::
getOptionListNext(std::string keyString, Value &value) const
{
  std::string pathString = this->getStringPath();
  PythonUtility::getOptionListNext<Value>(this->pythonConfig_, keyString, pathString, value);
}

//! extract a vector with unknown number of nEntries, must be a list
template<typename T>
void PythonConfig::
getOptionVector(std::string keyString, std::vector<T> &values) const
{
  std::string pathString = this->getStringPath();
  PythonUtility::getOptionVector(this->pythonConfig_, keyString, pathString, values);
}

//! extract a vector with known number of nEntries
template<typename T>
void PythonConfig::
getOptionVector(std::string keyString, int nEntries, std::vector<T> &values) const
{
  std::string pathString = this->getStringPath();

  PyObject *pyLocalValues = nullptr;
  pyLocalValues = getOptionPyObject(keyString, pyLocalValues);
  if (!pyLocalValues)
  {
    LOG(WARNING) << pathString << "[\"" << keyString << "\"]: no vector was given, set to zeros";
    values.resize(nEntries, T{0.0});
    return;
  }
  values = PythonUtility::convertFromPython<std::vector<T>>::get(pyLocalValues, T{0.0});
  if (values.size() < nEntries)
  {
    LOG(WARNING) << pathString << "[\"" << keyString << "\"]: given vector has only " << values.size() << " entries, fill with 0's to size "
      << nEntries;
    values.resize(nEntries, T{0.0});
  }
}
