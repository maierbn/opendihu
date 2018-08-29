#include "utility/string_utility.h"

#include <iostream>
#include <algorithm>
#include <iomanip>

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

//! replace from by to
std::string replace(std::string str, const std::string& from, const std::string& to)
{
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return str;
  std::string result(str);
  result.replace(start_pos, from.length(), to);
  return result;
}

template<>
std::string multiply<1>(std::string str)
{
  return str;
}

template<>
std::string multiply<2>(std::string str)
{
  return str+"*"+str;
}

template<>
std::string multiply<3>(std::string str)
{
  return str+"*"+str+"*"+str;
}

std::string extractBasename(std::string str)
{
  if (str.rfind(".") != std::string::npos)
  {
    str = str.substr(0, str.rfind("."));
  }
  if (str.find("/") != std::string::npos)
  {
    str = str.substr(str.rfind("/")+1);
  }
  return str;
}

}; // namespace
