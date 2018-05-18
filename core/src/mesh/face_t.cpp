#include "mesh/face_t.h"

#include "easylogging++.h"

namespace Mesh
{

face_t parseFace(std::string str)
{
  std::string faceNames[6] = {"0-", "0+", "1-", "1+", "2-", "2+"};

  for (int i = 0; i < 6; i++)
  {
    if (str == faceNames[i])
    {
      return face_t(i);
    }
  }

  LOG(ERROR) << "Could not parse face \"" << str << "\", possible values are: \"0-\", \"0+\", \"1-\", \"1+\", \"2-\", \"2+\"";
  return face_t::face0Minus;
}

template<>
Vec3 getNormal<3>(face_t face)
{
  switch(face)
  {
  case face0Minus:
    return Vec3({-1., 0., 0.});
  case face0Plus:
    return Vec3({1., 0., 0.});
  case face1Minus:
    return Vec3({0., -1, 0.});
  case face1Plus:
    return Vec3({0., 1., 0.});
  case face2Minus:
    return Vec3({0., 0., -1.});
  case face2Plus:
    return Vec3({0., 0., 1.});
  }
  return Vec3();
}

template<>
Vec2 getNormal<2>(face_t face)
{
  switch(face)
  {
  case face0Minus:
    return Vec2({-1., 0.});
  case face0Plus:
    return Vec2({1., 0.});
  case face1Minus:
    return Vec2({0., -1});
  case face1Plus:
    return Vec2({0., 1.});
  default:
    LOG(ERROR) << "Face not valid for 2D element in getNormal<2>";
  }
  return Vec2();
}

Vec3 getXiOnFace(face_t face, std::array<double,2> xiSurface)
{
  Vec3 xi;

  // xiSurface is 2D, coordinates on the face to integrate
  // set value of xi with 3D coordinates
  switch(face)
  {
  case face_t::face0Minus:
    xi = {0.0, xiSurface[1], xiSurface[0]};  // zy
    break;
  case face_t::face0Plus:
    xi = {1.0, xiSurface[0], xiSurface[1]};  // yz
    break;
  case face_t::face1Minus:
    xi = {xiSurface[0], 0.0, xiSurface[1]};  // xz
    break;
  case face_t::face1Plus:
    xi = {xiSurface[1], 1.0, xiSurface[0]};  // zx
    break;
  case face_t::face2Minus:
    xi = {xiSurface[1], xiSurface[0], 0.0};  // yx
    break;
  case face_t::face2Plus:
    xi = {xiSurface[0], xiSurface[1], 1.0};  // xy
    break;
  }
  return xi;
}

Vec2 getXiOnFace(face_t face, std::array<double,1> xiSurface)
{
  Vec2 xi;

  // xiSurface is 2D, coordinates on the face to integrate
  // set value of xi with 3D coordinates
  switch(face)
  {
  case face_t::face0Minus:
    xi = {0.0, xiSurface[0]};
    break;
  case face_t::face0Plus:
    xi = {1.0, xiSurface[0]};
    break;
  case face_t::face1Minus:
    xi = {xiSurface[0], 0.0};  // xz
    break;
  case face_t::face1Plus:
    xi = {xiSurface[0], 1.0};  // zx
    break;
  default:
    LOG(ERROR) << "Face not valid for 2D element in getXiOnFace";
  }
  return xi;
}

}  // namespace