#pragma once

#include <Python.h>  // has to be the first included header

#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <mutex>
#include <omp.h>

#include "control/types.h"

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
  static bool hasKey(const PyObject *settings, std::string key);

  //! checks if this settings is the empty list or None
  static bool isEmpty(const PyObject *settings, std::string key);

  //! checks if the object is a python list
  static bool isTypeList(const PyObject *object);

  //! given a python dictionary in settings, extract the value of given key and check if it is again a dict. Returns NULL, if the key does not exist. Then also a warning is printed.
  static PyObject *getOptionPyObject(const PyObject *settings, std::string key, std::string pathString, PyObject *defaultValue = Py_None);

  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  static double getOptionDouble(const PyObject *settings, std::string key, std::string pathString, double defaultValue, ValidityCriterion validityCriterion = None);

  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  static int getOptionInt(const PyObject *settings, std::string key, std::string pathString, int defaultValue, ValidityCriterion validityCriterion = None);

  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  static bool getOptionBool(const PyObject *settings, std::string key, std::string pathString, bool defaultValue = true);

  //! return the option value given by key in the python dictionary settings. If not found, return the defaultValue
  static std::string getOptionString(const PyObject *settings, std::string key, std::string pathString, std::string defaultValue);

  //! return the option value given by key in the python dictionary settings. If not found, return NULL
  static PyObject *getOptionFunction(const PyObject *settings, std::string key, std::string pathString);

  //! return the option value as array given by key in the python dictionary settings. If not found, return the defaultValue, also check if validityCriterion is met
  template<typename ValueType, int D>
  static std::array<ValueType, D> getOptionArray(PyObject* settings, std::string keyString, std::string pathString, std::array<ValueType, D> defaultValue,
                                                ValidityCriterion validityCriterion = None);

  //! return the option value as array given by key in the python dictionary settings. If not found, return the defaultValue, also check if validityCriterion is met
  template<typename ValueType, int D>
  static std::array<ValueType, D> getOptionArray(PyObject* settings, std::string keyString, std::string pathString, ValueType defaultValue,
                                                ValidityCriterion validityCriterion = None);

  //! Consider a Python dictionary in settings with key keyString. Set an internal iterator to the beginning and return the first key,value pair
  template<typename Key, typename Value>
  static std::pair<Key, Value> getOptionDictBegin(const PyObject *settings, std::string keyString, std::string pathString);

  //! If the internal iterator on the current dictionary is at the end of the dictionary
  static bool getOptionDictEnd(const PyObject *settings, std::string keyString, std::string pathString);

  //! Increment the internal iterator of the dictionary and set the next key,value pair in nextPair
  template<typename Key, typename Value>
  static void getOptionDictNext(const PyObject *settings, std::string keyString, std::string pathString, std::pair<Key, Value> &nextPair);

  //! Consider a Python list in settings with key keyString. Set an internal iterator to the beginning and return the first key,value pair
  template<typename Value>
  static Value getOptionListBegin(const PyObject *settings, std::string keyString, std::string pathString);

  //! If the internal iterator on the current list is at the end of the list
  static bool getOptionListEnd(const PyObject *settings, std::string keyString, std::string pathString);

  //! Increment the internal iterator of the list and set the next key,value pair in nextPair
  template<typename Value>
  static void getOptionListNext(const PyObject *settings, std::string keyString, std::string pathString, Value &value);

  //! extract a vector with exactly the specified number of nEntries, can be a dict or list, not specified entries are set to 0
  static void getOptionVector(const PyObject *settings, std::string keyString, std::string pathString, int nEntries, std::vector<double> &values);

  //! extract a vector with unknown number of nEntries, must be a list
  template<typename ValueType>
  static void getOptionVector(const PyObject *settings, std::string keyString, std::string pathString, std::vector<ValueType> &values);

  //! recursively print python dictionary to VLOG(1)
  static void printDict(PyObject *dict);

  //! recursively print a single python value
  static std::string getString(PyObject *value, int indent=0, int first_indent=0);

  template<typename T>
  struct convertFromPython
  {
    //! convert a python object to its corresponding c type, with type checking, if conversion is not possible, use defaultValue
    static T get(PyObject *object, T defaultValue);

    //! convert a python object to its corresponding c type, with type checking, if conversion is not possible use trivial default value (0 or 0.0 or "")
    static T get(PyObject *object);
  };

  template<typename T>
  struct convertToPython
  {
    //! convert a type to a python object
    static PyObject *get(T value);
  };

  //! create a python list out of the double vector
  static PyObject *convertToPythonList(std::vector<double> &data);

  //! create a python list out of the long vector
  static PyObject *convertToPythonList(std::vector<long> &data);

  //! create a python list out of the int vector
  static PyObject *convertToPythonList(std::vector<int> &data);

  //! create a python list out of the int vector
  static PyObject *convertToPythonList(std::vector<global_no_t> &data);

  //! create a python list out of the bool vector
  static PyObject *convertToPythonList(std::vector<bool> &data);

  //! create a python list out of the long array
  template<int D>
  static PyObject *convertToPythonList(std::array<long,D> &data);

  //! create a python list out of the bool array
  template<int D>
  static PyObject *convertToPythonList(std::array<bool,D> &data);

  //! create a python list out of the long vector
  static PyObject *convertToPythonList(unsigned int nEntries, double *data);

  //! create a python list from a double *
  static PyObject *convertToPythonList(double *value, int nValues);

  //! convert a PyUnicode object to a std::string
  static std::string pyUnicodeToString(PyObject *object);

  //! print a error message to stderr if there was a python error
  static void checkForError();

  /** Helper class that acquires the global interpreter lock of the python interpreter. This is needed for every call to the python api when running multi-threaded (openmp) programs.
   * A critical section for python API code starts when an object of this class is instantiated and ends, when the object gets destructed.
   * Typical usage :
   * {
   *   GlobalInterpreterLock lock;
   *   // python c-api code here
   * }
   */
  class GlobalInterpreterLock
  {
  public:
    //! constructor
    GlobalInterpreterLock();
    
    //! destructor
    ~GlobalInterpreterLock();
  private:
    PyGILState_STATE gstate_;
    //static int nGILS_;
    //static std::map<int, int> nGilsThreads_;
    //PyThreadState *mainThreadState_;
    
    //static std::recursive_mutex mutex_;  //< mutex for critical section
    //static std::unique_lock<std::recursive_mutex> lock_;
    
    //static bool lockInitialized_;
    //static omp_nest_lock_t lock_;
  };
  
private:

  static PyObject *itemList;    //< list of items (key,value) for dictionary,  to use for getOptionDictBegin, getOptionDictEnd, getOptionDictNext
  static int itemListIndex;     //< current index of itemList

  static PyObject *list;        //< python list to use for getOptionListBegin, getOptionListEnd, getOptionListNext
  static int listIndex;         //< current index for list
};

//! output python object
std::ostream &operator<<(std::ostream &stream, PyObject *object);


#include "utility/python_utility.tpp"
#include "utility/python_utility_convert_scalar.tpp"
#include "utility/python_utility_convert_lists.tpp"
