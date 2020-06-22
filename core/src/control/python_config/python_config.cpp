#include "control/python_config/python_config.h"

#include <sstream>

#include "easylogging++.h"
#include "utility/python_utility.h"

//! constructor from python object
PythonConfig::PythonConfig(PyObject *specificSettings)
{
  pythonConfig_ = specificSettings;
  VLOG(1) << "PythonConfig::constructor " << PythonUtility::getString(pythonConfig_);
  Py_XINCREF(pythonConfig_);
}

//! constructor as sub scope of another python config
PythonConfig::PythonConfig(const PythonConfig &rhs, std::string key)
{
  pythonConfig_ = rhs.getOptionPyObject(key);
  VLOG(1) << "PythonConfig::constructor(rhs,key=\"" << key << "\") " << PythonUtility::getString(pythonConfig_);
  Py_XINCREF(pythonConfig_);

  int pathSize = std::distance(rhs.pathBegin(), rhs.pathEnd());
  path_.resize(pathSize+1);
  std::copy(rhs.pathBegin(), rhs.pathEnd(), path_.begin());
  path_[pathSize] = key;
}

//! constructor as sub scope of another python config which is a list
PythonConfig::PythonConfig(const PythonConfig &rhs, int i)
{
  // store updated path
  int pathSize = std::distance(rhs.pathBegin(), rhs.pathEnd());
  path_.resize(pathSize+1);
  std::copy(rhs.pathBegin(), rhs.pathEnd(), path_.begin());
  std::stringstream s;
  s << i;
  path_[pathSize] = std::string(s.str());

  // check if the current object is a list
  if (PyList_Check(rhs.pyObject()))
  {
    int nEntries = PyList_Size(rhs.pyObject());
    if (i >= nEntries)
    {
      LOG(ERROR) << getStringPath() << " list has only " << nEntries << " entries, but entry "
        << i << " is required.";
      pythonConfig_ = PyList_GetItem(rhs.pyObject(), (Py_ssize_t)nEntries-1);
    }
    else
    {
      pythonConfig_ = PyList_GetItem(rhs.pyObject(), (Py_ssize_t)i);
    }
  }
  else
  {
    LOG(WARNING) << getStringPath() << " is not a list";
    pythonConfig_ = rhs.pyObject();
  }

}

//! constructor directly from PyObject*, path from rhs + key
PythonConfig::PythonConfig(const PythonConfig &rhs, std::string key, PyObject *config)
{
  pythonConfig_ = config;
  VLOG(1) << "PythonConfig::constructor(rhs,key=\"" << key << "\",config) " << PythonUtility::getString(pythonConfig_);
  Py_XINCREF(pythonConfig_);

  int pathSize = std::distance(rhs.pathBegin(), rhs.pathEnd());
  path_.resize(pathSize+2);
  std::copy(rhs.pathBegin(), rhs.pathEnd(), path_.begin());
  path_[pathSize] = key;
  path_[pathSize+1] = std::string("...");
}

//! constructor directly from PyObject*, path from rhs + key
PythonConfig::PythonConfig(const PythonConfig &rhs, std::string key, std::string key2, PyObject *config)
{
  pythonConfig_ = config;
  VLOG(1) << "PythonConfig::constructor(rhs,key=\"" << key << "\",config) " << PythonUtility::getString(pythonConfig_);
  Py_XINCREF(pythonConfig_);

  int pathSize = std::distance(rhs.pathBegin(), rhs.pathEnd());
  path_.resize(pathSize+2);
  std::copy(rhs.pathBegin(), rhs.pathEnd(), path_.begin());
  path_[pathSize] = key;
  path_[pathSize+1] = key2;
}

PythonConfig::~PythonConfig()
{
  //Py_XDECREF(pythonConfig_);
}

PythonConfig &PythonConfig::operator=(const PythonConfig &rhs)
{
  pythonConfig_ = rhs.pythonConfig_;
  VLOG(1) << "PythonConfig::operator=(rhs)" << PythonUtility::getString(pythonConfig_);
  Py_XINCREF(pythonConfig_);

  int pathSize = std::distance(rhs.pathBegin(), rhs.pathEnd());
  path_.resize(pathSize);
  std::copy(rhs.pathBegin(), rhs.pathEnd(), path_.begin());
  return *this;
}

//! return the actual PyObject object
PyObject *PythonConfig::pyObject() const
{
  return pythonConfig_;
}

void PythonConfig::setPyObject(PyObject *pyObject)
{
  pythonConfig_ = pyObject;
  Py_XINCREF(pythonConfig_);
}

std::vector<std::string>::const_iterator PythonConfig::pathBegin() const
{
  return path_.begin();
}

std::vector<std::string>::const_iterator PythonConfig::pathEnd() const
{
  return path_.end();
}

std::string PythonConfig::getStringPath() const
{
  std::stringstream pathString;
  pathString << "config";
  for (std::vector<std::string>::const_iterator iter = path_.begin(); iter != path_.end(); iter++)
  {
    if (*iter == "...")
    {
      pathString << "[...]";
    }
    else
    {
      pathString << "[\"" << *iter << "\"]";
    }
  }
  return pathString.str();
}

//! checks if the settings contain the given key, no warning is printed
bool PythonConfig::hasKey(std::string key) const
{
  return PythonUtility::hasKey(this->pythonConfig_, key);
}

//! checks if the settings contain the given key, no warning is printed
bool PythonConfig::isEmpty(std::string key) const
{
  return PythonUtility::isEmpty(this->pythonConfig_, key);
}

//! return all keys of the current dict as vector of strings
void PythonConfig::getKeys(std::vector<std::string> &keys)
{
  if (!PyDict_Check(pythonConfig_))
    return;

  PyObject *pyList = PyDict_Items(pythonConfig_);
  std::vector<std::pair<std::string,PyObject *>> list = PythonUtility::convertFromPython<std::vector<std::pair<std::string,PyObject *>>>::get(pyList);

  keys.clear();
  for (std::vector<std::pair<std::string,PyObject *>>::iterator iter = list.begin(); iter != list.end(); iter++)
  {
    keys.push_back(iter->first);
  }
}

//! given a python dictionary in settings, extract the value of given key and check if it is again a dict. Returns NULL, if the key does not exist. Then also a warning is printed.
PyObject *PythonConfig::getOptionPyObject(std::string key, PyObject *defaultValue) const
{
  std::string pathString = getStringPath();
  return PythonUtility::getOptionPyObject(this->pythonConfig_, key, pathString, defaultValue);
}

//! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
double PythonConfig::getOptionDouble(std::string key, double defaultValue, PythonUtility::ValidityCriterion validityCriterion) const
{
  std::string pathString = getStringPath();
  return PythonUtility::getOptionDouble(this->pythonConfig_, key, pathString, defaultValue, validityCriterion);
}

//! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
int PythonConfig::getOptionInt(std::string key, int defaultValue, PythonUtility::ValidityCriterion validityCriterion) const
{
  std::string pathString = getStringPath();
  return PythonUtility::getOptionInt(this->pythonConfig_, key, pathString, defaultValue, validityCriterion);
}

//! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
bool PythonConfig::getOptionBool(std::string key, bool defaultValue) const
{
  std::string pathString = getStringPath();
  return PythonUtility::getOptionBool(this->pythonConfig_, key, pathString, defaultValue);
}

//! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
std::string PythonConfig::getOptionString(std::string key, std::string defaultValue) const
{
  std::string pathString = getStringPath();
  return PythonUtility::getOptionString(this->pythonConfig_, key, pathString, defaultValue);
}

//! return the option value given by key in the python dictionary settings. If not found, return NULL
PyObject *PythonConfig::getOptionFunction(std::string key) const
{
  std::string pathString = getStringPath();
  return PythonUtility::getOptionFunction(this->pythonConfig_, key, pathString);
}

//! If the internal iterator on the current dictionary is at the end of the dictionary
bool PythonConfig::getOptionDictEnd(std::string keyString) const
{
  std::string pathString = getStringPath();
  return PythonUtility::getOptionDictEnd(this->pythonConfig_, keyString, pathString);
}

//! If the internal iterator on the current list is at the end of the list
bool PythonConfig::getOptionListEnd(std::string keyString) const
{
  std::string pathString = getStringPath();
  return PythonUtility::getOptionListEnd(this->pythonConfig_, keyString, pathString);
}

//! extract a vector with exactly the specified number of nEntries, can be a dict or list, not specified entries are set to 0
void PythonConfig::getOptionVector(std::string keyString, int nEntries, std::vector<double> &values) const
{
  std::string pathString = getStringPath();
  PythonUtility::getOptionVector(this->pythonConfig_, keyString, pathString, nEntries, values);
}

void PythonConfig::
getOptionVector(std::string keyString, std::vector<PythonConfig> &values) const
{
  PyObject *settings = this->pythonConfig_;
  if (settings)
  {
    // check if input dictionary contains the key
    PyObject *key = PyUnicode_FromString(keyString.c_str());
    if (PyDict_Contains((PyObject *)settings, key))
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
        PyObject *item = this->getOptionListBegin<PyObject *>(keyString);

        // loop over other values
        for (;
            !this->getOptionListEnd(keyString);
            this->getOptionListNext<PyObject *>(keyString, item))
        {
          PythonConfig configItem(*this, keyString, item);
          values.push_back(configItem);
        }
      }
      else
      {
        // not a list, but a single entry (only 1 entry)
        PyObject *item = (PyObject *)settings;
        PythonConfig configItem(*this, keyString, item);
        values.push_back(configItem);
      }
    }
    else
    {
      // this is no warning
      std::string pathString = this->getStringPath();
      LOG(DEBUG) << "" << pathString << "[\"" << keyString << "\"] not set in \"" << Control::settingsFileName << "\".";
    }
    Py_CLEAR(key);
  }
}

std::ostream &operator<<(std::ostream &stream, const PythonConfig rhs)
{
  if (rhs.pyObject() == nullptr)
  {
    stream << std::string("(null)");
  }
  else
  {
    stream << rhs.getStringPath();
  }
  return stream;
}
