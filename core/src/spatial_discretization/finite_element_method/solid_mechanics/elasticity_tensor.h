#pragma once

#include <Python.h>  // has to be the first included header

#include <array>

namespace SpatialDiscretization
{
 
/** Class to store an elasticity tensor for 3D continuum mechanics.
 * Due to symmetrics C_ABCD = C_BACD = C_ABDC, is has 21 distinct entries.
 */
class ElasticityTensor : 
  public std::array<double, 21>
{
 
public:
 
  //! return the entry klrs of the elasticity tensor
  double getEntry(int k, int l, int r, int s);
  
protected:
 
  //! return the index to the elastcity array with 21 distinct entries, for logical component klrs
  int getEntryNo(int k, int l, int r, int s);
 
};

};  // namespace

