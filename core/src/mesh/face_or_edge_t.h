#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>

#include "control/types.h"

namespace Mesh
{

/** An enum that is compatible to face_t and adds four edges (thought in 3D).
 *  These values are used to index ghost meshes around a domain.
 */
enum face_or_edge_t
{
  faceEdge0Minus = 0, faceEdge0Plus,
  faceEdge1Minus, faceEdge1Plus,
  faceEdge2Minus, faceEdge2Plus,
  edge0Minus1Minus, edge0Plus1Minus,
  edge0Minus1Plus,  edge0Plus1Plus
};

//! get the opposite face, i.e. edge0Minus1Minus <--> edge0Plus1Plus, edge0Plus1Minus <--> edge0Minus1Plus
face_or_edge_t oppositeEdgeOrFace(face_or_edge_t face);

//! get the string "0-", "0+", "1-", "1+", "2-", "2+", "0-1-", "0+1-", "0-1+", "0+1+"
std::string getString(face_or_edge_t face);

}  // namespace

