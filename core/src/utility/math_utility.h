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
  static double distance(const Vec3 node1, const Vec3 node2);
  
  //! return the euclidean norm of the vector
  static double length(const Vec3 node);
  
  //! same as length, return the euclidean norm of the vector
  static double norm(const Vec3 node);
  
  //! compute the matrix T = J^{-1}J^{-T} when J is the jacobian
  //! returns the matrix in row major storage order, matrix is symmetric and the determinant of the jacobian
  static std::array<double,9> computeTransformationMatrixAndDeterminant(const std::array<Vec3,3> &jacobian, double &determinant);
  
  //! computes the determinant of the matrix 
  static double computeDeterminant(const std::array<Vec3,3> &jacobian);
  
  //! computes v1^T * T * v2 where T is the symmetric transformation matrix (3D)
  static double applyTransformation(const std::array<double,9> &transformationMatrix, const Vec3 &vector1, const Vec3 &vector2);
  
  //! computes v1^T * T * v2 where T is the symmetric transformation matrix (2D)
  static double applyTransformation(const std::array<double,4> &transformationMatrix, const Vec2 &vector1, const Vec2 &vector2);  
  
  //! computes the factor J_D that is needed when transforming an integral from world coordinates to parameter space
  template<int D>
  static double computeIntegrationFactor(const std::array<Vec3,D> &jacobian);
  
  static constexpr auto PI = 3.14159265358979323846;
  
  //! check if vector b is a subsequence of a, i.e. is contained in a
  static bool isSubsequenceOf(std::vector<int> a, std::vector<int> b, int &subsequenceAStartPos);
};

Vec3 operator-(Vec3 node1, Vec3 node2);
Vec3 operator+(Vec3 node1, Vec3 node2);
Vec3 &operator+=(Vec3 &node1, Vec3 node2);
Vec3 operator*(double lambda, Vec3 node);
template<unsigned long N>
std::ostream &operator<<(std::ostream &stream, const std::array<double,N> node);
template<typename T>
bool operator==(const std::vector<T> &vec1, const std::vector<T> &vec2);

#include "utility/math_utility.tpp"