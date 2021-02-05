
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <sstream>

int main(int argc, char *argv[])
{

  std::stringstream compileCommand;
  compileCommand << "gcc src/hodgkin_huxley_1952_1_0___gpu_101.0.c -fPIC -O3 -march=native -shared-fopenmp-foffload=\"-O3 -lm blabla\"-o lib/hodgkin_huxley_1952_1_0___gpu_101.so.0  && mv lib/hodgkin_huxley_1952_1_0___gpu_101.so.0 lib/hodgkin_huxley_1952_1_0___gpu_101.so";

  std::cout << "compileCommand: \"" << compileCommand.str() << "\"." << std::endl;
  
  // remove "-fopenmp" in the compile command
  std::string newCompileCommand = compileCommand.str();
  std::string strToReplace = "-fopenmp";
  std::size_t pos = newCompileCommand.find(strToReplace);
  newCompileCommand.replace(pos, strToReplace.length(), "");

  std::cout << "newCompileCommand: \"" << newCompileCommand << "\"." << std::endl;
  
  // remove -foffload="..."strToReplace = "-fopenmp";
  pos = newCompileCommand.find("-foffload=\"");
  std::size_t pos2 = newCompileCommand.find("\"", pos+11);
  std::cout << "pos=" << pos << ",  pos2=" << pos2 << std::endl;
  newCompileCommand.replace(pos, pos2-pos+1, "");

  std::cout << "newCompileCommand: \"" << newCompileCommand << "\"." << std::endl;
  
  return EXIT_SUCCESS;
}
