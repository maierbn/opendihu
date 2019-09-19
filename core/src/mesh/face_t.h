#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>

#include "control/types.h"

namespace Mesh
{

enum face_t
{
  face0Minus = 0, face0Plus,
  face1Minus, face1Plus,
  face2Minus, face2Plus
};

//! parse the face value from a string
face_t parseFace(std::string str);

//! get the opposite face, i.e. face0Minus <--> face0Plus, face1Minus <--> face1Plus, face2Minus <--> face2Plus
face_t oppositeFace(face_t face);

//! get the string "0-", "0+", "1-", "1+", "2-" or "2+"
std::string getString(face_t face);

//! return the outward normal vector to the face
template<int D>
std::array<double,D> getNormal(face_t face);

//! return a value of xi in [0,1]^3 that lies on the face in parameter space, the position on the face is given by xiSurface in [0,1]^2
Vec3 getXiOnFace(face_t face, std::array<double,2> xiSurface);

//! return a value of xi in [0,1]^2 that lies on the face in parameter space, the position on the face is given by xiSurface in [0,1]
Vec2 getXiOnFace(face_t face, std::array<double,1> xiSurface);

//! dummy function
VecD<1> getXiOnFace(face_t face, std::array<double,0> xiSurface);

}  // namespace

