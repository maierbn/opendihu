#include "spatial_discretization/finite_element_method/solid_mechanics/elasticity_tensor.h"

#include <Python.h>  // has to be the first included header

#include "easylogging++.h"

namespace SpatialDiscretization
{

//! return the entry klrs of the elasticity tensor
int ElasticityTensor::
getEntryNo(int k, int l, int r, int s)
{
  // this method was tested outside of this codebase and is correct
  // the inverse mapping is given by the following array:
  //int indices[21][4] = {
  //  {0,0,0,0},{0,1,0,0},{0,2,0,0},{1,1,0,0},{1,2,0,0},{2,2,0,0},{0,1,0,1},{0,2,0,1},{1,1,0,1},{1,2,0,1},
  //  {2,2,0,1},{0,2,0,2},{1,1,0,2},{1,2,0,2},{2,2,0,2},{1,1,1,1},{1,2,1,1},{2,2,1,1},{1,2,1,2},{2,2,1,2},
  //  {2,2,2,2}
  //};
  switch(k)
  {
  case 0: // k
    switch(l)
    {
    case 0: // l
      switch(r)
      {
      case 0: // 000
        switch(s)
        {
        case 0:
          return 0;
          break;
        case 1:
          return 1;
          break;
        case 2:
          return 2;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 001
        switch(s)
        {
        case 0:
          return 1;
          break;
        case 1:
          return 3;
          break;
        case 2:
          return 4;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 002
        switch(s)
        {
        case 0:
          return 2;
          break;
        case 1:
          return 4;
          break;
        case 2:
          return 5;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 1: // l
      switch(r)
      {
      case 0: // 010
        switch(s)
        {
        case 0:
          return 1;
          break;
        case 1:
          return 6;
          break;
        case 2:
          return 7;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 011
        switch(s)
        {
        case 0:
          return 6;
          break;
        case 1:
          return 8;
          break;
        case 2:
          return 9;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 012
        switch(s)
        {
        case 0:
          return 7;
          break;
        case 1:
          return 9;
          break;
        case 2:
          return 10;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 2: // l
      switch(r)
      {
      case 0: // 020
        switch(s)
        {
        case 0:
          return 2;
          break;
        case 1:
          return 7;
          break;
        case 2:
          return 11;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 021
        switch(s)
        {
        case 0:
          return 7;
          break;
        case 1:
          return 12;
          break;
        case 2:
          return 13;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 022
        switch(s)
        {
        case 0:
          return 11;
          break;
        case 1:
          return 13;
          break;
        case 2:
          return 14;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    default:
      LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
      break;
    }  // l
    break;
  case 1: // k
    switch(l)
    {
    case 0: // l
      switch(r)
      {
      case 0: // 100
        switch(s)
        {
        case 0:
          return 1;
          break;
        case 1:
          return 6;
          break;
        case 2:
          return 7;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 101
        switch(s)
        {
        case 0:
          return 6;
          break;
        case 1:
          return 8;
          break;
        case 2:
          return 9;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 102
        switch(s)
        {
        case 0:
          return 7;
          break;
        case 1:
          return 9;
          break;
        case 2:
          return 10;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 1: // l
      switch(r)
      {
      case 0: // 110
        switch(s)
        {
        case 0:
          return 3;
          break;
        case 1:
          return 8;
          break;
        case 2:
          return 12;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 111
        switch(s)
        {
        case 0:
          return 8;
          break;
        case 1:
          return 15;
          break;
        case 2:
          return 16;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 112
        switch(s)
        {
        case 0:
          return 12;
          break;
        case 1:
          return 16;
          break;
        case 2:
          return 17;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 2: // l
      switch(r)
      {
      case 0: // 120
        switch(s)
        {
        case 0:
          return 4;
          break;
        case 1:
          return 9;
          break;
        case 2:
          return 13;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 121
        switch(s)
        {
        case 0:
          return 9;
          break;
        case 1:
          return 16;
          break;
        case 2:
          return 18;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 122
        switch(s)
        {
        case 0:
          return 13;
          break;
        case 1:
          return 18;
          break;
        case 2:
          return 19;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    default:
      LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
      break;
    }  // l
    break;
  case 2: // k
    switch(l)
    {
    case 0: // l
      switch(r)
      {
      case 0: // 200
        switch(s)
        {
        case 0:
          return 2;
          break;
        case 1:
          return 7;
          break;
        case 2:
          return 11;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 201
        switch(s)
        {
        case 0:
          return 7;
          break;
        case 1:
          return 12;
          break;
        case 2:
          return 13;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 202
        switch(s)
        {
        case 0:
          return 11;
          break;
        case 1:
          return 13;
          break;
        case 2:
          return 14;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 1: // l
      switch(r)
      {
      case 0: // 210
        switch(s)
        {
        case 0:
          return 4;
          break;
        case 1:
          return 9;
          break;
        case 2:
          return 13;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 211
        switch(s)
        {
        case 0:
          return 9;
          break;
        case 1:
          return 16;
          break;
        case 2:
          return 18;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 212
        switch(s)
        {
        case 0:
          return 13;
          break;
        case 1:
          return 18;
          break;
        case 2:
          return 19;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 2: // l
      switch(r)
      {
      case 0: // 220
        switch(s)
        {
        case 0:
          return 5;
          break;
        case 1:
          return 10;
          break;
        case 2:
          return 14;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 221
        switch(s)
        {
        case 0:
          return 10;
          break;
        case 1:
          return 17;
          break;
        case 2:
          return 19;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 222
        switch(s)
        {
        case 0:
          return 14;
          break;
        case 1:
          return 19;
          break;
        case 2:
          return 20;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    default:
      LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
      break;
    }  // l
    break;
  default:
    LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
    break;
  }  // k
  return 0;
}

//! return the entry klrs of the elasticity tensor
double ElasticityTensor::
getEntry(int k, int l, int r, int s)
{
  return this->operator[](getEntryNo(k,l,r,s));
}

};    // namespace