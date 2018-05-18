#include "utility/string_utility.h"

#include <iostream>
#include <algorithm>
#include <iomanip>

namespace StringUtility
{

template<typename IterType>
void outputValuesBlock(std::ostream &stream, IterType valuesBegin,
                       IterType valuesEnd, int nValuesPerRow)
{
  int i = 0;
  std::size_t nValues = valuesEnd - valuesBegin;
  for (IterType valuesIter = valuesBegin; valuesIter != valuesEnd; valuesIter++, i++)
  {
    stream << (*valuesIter >= 0? " " : "") << " " << std::scientific << std::setw(7) << *valuesIter;

    // add newline after nValuesPerRow entries per row
    if (nValuesPerRow != -1 && (i+1) % nValuesPerRow == 0 && i < nValues-1)
    {
      stream << std::endl;
    }
  }
  stream << std::endl;
};


}; // namespace
