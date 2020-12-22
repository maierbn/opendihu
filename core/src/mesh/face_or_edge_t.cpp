#include "mesh/face_or_edge_t.h"

#include "easylogging++.h"

namespace Mesh
{

std::string getString(face_or_edge_t face)
{
  switch(face)
  {
  case faceEdge0Minus:
    return std::string("0-");
  case faceEdge0Plus:
    return std::string("0+");
  case faceEdge1Minus:
    return std::string("1-");
  case faceEdge1Plus:
    return std::string("1+");
  case faceEdge2Minus:
    return std::string("2-");
  case faceEdge2Plus:
    return std::string("2+");
  case edge0Minus1Minus:
    return std::string("0-1-");
  case edge0Plus1Minus:
    return std::string("0+1-");
  case edge0Minus1Plus:
    return std::string("0-1+");
  case edge0Plus1Plus:
    return std::string("0+1+");
  }
  return std::string("");
}
face_or_edge_t oppositeFace(face_or_edge_t face)
{
  switch(face)
  {
  case faceEdge0Minus:
  case faceEdge1Minus:
  case faceEdge2Minus:
    return face_or_edge_t((int)face + 1);
  case faceEdge0Plus:
  case faceEdge1Plus:
  case faceEdge2Plus:
    return face_or_edge_t((int)face - 1);
  case edge0Minus1Minus:
    return edge0Plus1Plus;
  case edge0Plus1Minus:
    return edge0Minus1Plus;
  case edge0Minus1Plus:
    return edge0Plus1Minus;
  case edge0Plus1Plus:
    return edge0Minus1Minus;
  }
#ifndef __PGI
  return faceEdge0Minus;
#endif
}

}  // namespace
