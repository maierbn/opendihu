#pragma once

#include <petscsys.h>
#include <Python.h>
#include <list>
#include <memory>

#include <output_writer/generic.h>
#include <data_management/data.h>

class MeshManager;
class DihuContext
{
public:
  ///! constructor, initialize context, parse command line parameters and input file
  DihuContext(int argc, char *argv[]);
  
  ///! constructor for test cases
  DihuContext(int argc, char *argv[], std::string pythonSettings);
  
  ///! default copy-constructor
  DihuContext(const DihuContext &rhs) = default;
  DihuContext(DihuContext &&rhs) = default;
  
  ///! return a context object with config originated at child node with given key
  DihuContext operator[](std::string keyString) const;
  
  ///! return the top-level python config object
  PyObject *getPythonConfig() const;
  
  ///! return the mesh manager object that contains all meshes
  std::shared_ptr<MeshManager> meshManager() const;
  
  ///! call all output writers to write output, timeStepNo of -1 means no time step number in output filename
  void writeOutput(Data::Data &problemData, int timeStepNo = -1, double currentTime = 0.0) const;
  
  ///! destructor
  ~DihuContext();
  
private:
  ///! read in file and execute python script and store global variables
  void loadPythonScriptFromFile(std::string filename);
  
  ///! execute python script and store global variables
  void loadPythonScript(std::string text);
  
  ///! initiaize the easylogging++ framework
  void initializeLogging(int argc, char *argv[]);
  
  ///! create output writer from settings
  void initializeOutputWriter();
  
  ///! helper function that creates on outputWriter
  void createOutputWriterFromSettings(PyObject *dict);
  
  PyObject *pythonConfig_;    ///< the top level python config dictionary
  
  static std::shared_ptr<MeshManager> meshManager_;   ///< object that saves all meshes that are used
  static std::list<std::unique_ptr<OutputWriter::Generic>> outputWriter_;    ///< list of active output writers
  static bool initialized_;
};