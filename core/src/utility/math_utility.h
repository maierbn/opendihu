#pragma once

#include <vector>
#include <streambuf>

#include <petscmat.h>
#include "control/types.h"


class MathUtility
{
public:
  
  //! returns v*v (square)
  static double sqr(double v);
  
  //! returns v*v (square)
  static int sqr(int v);
  
  //! return the euclidean norm of the difference vector, i.e. the distance between the nodes
  static double distance(Vec3 node1, Vec3 node2);
  
  //! return the euclidean norm of the vector
  static double length(Vec3 node);
  
  //! compute the matrix T = J^{-1}J^{-T} when J is given as (jacobianColumn0 jacobianColumn1 jacobianColumn2)
  //! returns the matrix in row major storage order, matrix is symmetric and the determinant of the jacobian
  static std::array<double,9> computeTransformationMatrixAndDeterminant(
    Vec3 &jacobianColumn0, Vec3 &jacobianColumn1, Vec3 &jacobianColumn2, double &determinant);
  
  //! computes the determinant of the matrix 
  static double computeDeterminant(Vec3 &jacobianColumn0, Vec3 &jacobianColumn1, Vec3 &jacobianColumn2);
  
  //! computes v1^T * T * v2 where T is the symmetric transformation matrix (3D)
  static double applyTransformation(std::array<double,9> &transformationMatrix, Vec3 &vector1, Vec3 &vector2);
  
  //! computes v1^T * T * v2 where T is the symmetric transformation matrix (2D)
  static double applyTransformation(std::array<double,4> &transformationMatrix, Vec2 &vector1, Vec2 &vector2);
  
  //! computes the factor J_2 that is needed when transforming an integral from world coordinates to 2D parameter space
  static double compute2DIntegrationFactor(Vec3 &jacobianColumn0, Vec3 &jacobianColumn1);
  
  static constexpr auto PI = 3.14159265358979323846;
};

Vec3 operator-(Vec3 node1, Vec3 node2);
Vec3 operator+(Vec3 node1, Vec3 node2);
Vec3 &operator+=(Vec3 &node1, Vec3 node2);
Vec3 operator*(double lambda, Vec3 node);
std::ostream &operator<<(std::ostream &stream, Vec2 node);
std::ostream &operator<<(std::ostream &stream, Vec3 node);