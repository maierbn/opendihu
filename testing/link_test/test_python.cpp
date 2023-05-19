
#include <Python.h>
#include <iostream>

int main(int argc, char** argv)
{
  /* Code to test Python.h */
  Py_Initialize();
  std::cout << "Python initialized." << std::endl;
  Py_Finalize();
  std::cout << "Python finalized." << std::endl;
  
  return EXIT_SUCCESS;
}