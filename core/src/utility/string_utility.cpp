#include "utility/string_utility.h"

#include <iostream>
#include <algorithm>
#include <iomanip>
#ifdef __GNUC__
#include <cxxabi.h>
#endif

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
  str.erase(str.begin(), std::find_if (str.begin(), str.end(), [](char c){return !std::isspace(c);}));

  // remove whitespace at the end
  str.erase(std::find_if (str.rbegin(), str.rend(),
             [](char c){ return !std::isspace(c); }).base(),
            str.end());
}

//! replace from by to
std::string replace(std::string str, const std::string& from, const std::string& to)
{
  size_t start_pos = str.find(from);
  if (start_pos == std::string::npos)
    return str;
  std::string result(str);
  result.replace(start_pos, from.length(), to);
  return result;
}

//! replace from by to
std::string replaceAll(std::string str, const std::string& from, const std::string& to)
{
  std::string result(str);
  while(result.find(from) != std::string::npos)
  {
    result = replace(result, from, to);
  }
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

std::string timeToString(const tm* const time)
{   
  // to format: %Y/%m/%d %H:%M:%S
  std::string date;
    
  date += std::to_string( time->tm_year + 1900 ) + "/"
       +  std::to_string( time->tm_mon + 1 ) + "/"
       +  std::to_string( time->tm_mday ) + " ";
  if( time->tm_hour < 10 )
  {
      date += "0";
  }
  date += std::to_string( time->tm_hour ) + ":";
   if( time->tm_min < 10 )
  {
      date += "0";
  }
  date += std::to_string( time->tm_min ) + ":";
  if( time->tm_sec < 10 )
  {
      date += "0";
  }
  date += std::to_string( time->tm_sec );
  return date;
}

std::string demangle(const char *typeidName)
{
#ifdef __GNUC__
  // source: https://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html
  int status;
  return std::string(abi::__cxa_demangle(typeidName, 0, 0, &status));
#else
  return std::string(typeidName);
#endif
}

std::size_t stringLength(std::string string)
{
  int length = string.length();

  // unicode characters start with -62 and use 2 chars
  // or with -30 and use 3 chars
  length -= std::count(string.begin(), string.end(), char(-50));
  length -= std::count(string.begin(), string.end(), char(-62));
  length -= 2*std::count(string.begin(), string.end(), char(-30));

  return length;
}

}  // namespace
