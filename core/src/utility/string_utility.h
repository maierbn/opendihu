#pragma once

#include <string>
#include <vector>

namespace StringUtility
{

//! extract from a line in a file the literal that follows after a string key and convert it to int
int getNumberAfterString(std::string line, std::string key);

//! remove beginning of line until (including) key, return removed substring
std::string extractUntil(std::string &line, std::string key);

//! remove whitespace (' ', '\t', '\n') at the beginning and end of the string
void trim(std::string &str);

//! output the values separated by spaces, after nValuesPerRow there will be a line break
void outputValuesBlock(std::ostream &stream, const std::vector<double> &values, int nValuesPerRow);

};