#include "utility/string_utility.h"

#include <iostream>
#include <algorithm>

namespace StringUtility
{

int getNumberAfterString(std::string line, std::string key)
{
  std::string value = line.substr(line.find(key)+key.length());
  
  // remove leading white space and '=' signs
  while(isspace(value[0]) || value[0] == '=')
    value = value.substr(1);
  return atoi(value.c_str());
}

std::string extractUntil(std::string &line, std::string key)
{
  std::string result;
  if (line.find(key) != std::string::npos)
  {
    result = line.substr(0, line.find(key));
    line = line.substr(line.find(key) + key.length());
  }
  return result;
}

void trim(std::string &str)
{
  // remove whitespace at the beginning
  str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](char c){return !std::isspace(c);}));
  
  // remove whitespace at the end
  str.erase(std::find_if(str.rbegin(), str.rend(),
             [](char c){ return !std::isspace(c); }).base(),
            str.end());
}

void outputValuesBlock(std::ostream &stream, const std::vector<double> &values, int nValuesPerRow)
{
  for (unsigned int i = 0; i < values.size(); i++)
  {
    stream << " " << values[i];
    
    // add newline after nValuesPerRow entries per row
    if (i+1 % nValuesPerRow == 0)
    {
      stream << std::endl << " "; 
    }
  }
  stream << std::endl;
};

}; // namespace