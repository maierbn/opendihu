#include "utility/vector_operators.h"

#include <sstream>

#include "utility/petsc_utility.h"
#include "easylogging++.h"

template<>
std::ostream &operator<<(std::ostream &stream, const std::vector<double> &values)
{
  if (values.empty())
  {
    stream << "[]";
    return stream;
  }

  stream << "[";

  // first entry
  if (values[0] == std::numeric_limits<double>::max())
    stream << "None";
  else
    stream << values[0];

  if (VLOG_IS_ON(1))
  {
    // with VLOG output all entries
    for (unsigned long i = 1; i < values.size(); i++)
    {
      stream << ",";
      if (values[i] == std::numeric_limits<double>::max())
        stream << "None";
      else
        stream << values[i];
    }
  }
  else
  {
    // without VLOG only output the first 100 entries
    unsigned long i = 1;
    for (; i < std::min(100ul,values.size()); i++)
    {
      stream << ",";
      if (values[i] == std::numeric_limits<double>::max())
        stream << "None";
      else
        stream << values[i];
    }
    if (i == 100 && i < values.size())
      stream << "... " << values.size() << " entries total, only showing the first 100 (call with -vmodule=vector_operators*=1 to show all)";
  }

  stream << "]";
  return stream;
}
