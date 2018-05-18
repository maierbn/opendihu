#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"
#include "equation/mooney_rivlin_compressible.h"

namespace SpatialDiscretization
{


/** specialisation for compressible mooney rivlin
 * This is not yet fully implemented.
 */
/*
template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
class FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType, MixedQuadratureType, Equation::isCompressible<Term>, Mesh::isDeformable<typename BasisOnMeshType::Mesh>
> :
  public FiniteElementMethodBase<BasisOnMeshType, MixedQuadratureType>,
  public SolidMechanicsUtility<BasisOnMeshType, MixedQuadratureType, Equation::Static::CompressibleMooneyRivlin>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, MixedQuadratureType>::FiniteElementMethodBase;

  //! assemble the stiffness matrix
  void setStiffnessMatrix();

  //! initialize material constants from config
  void initialize() override;

private:

  //! assemble right hand side
  void manipulateWeakRhs() override;

  //! compute the reduced invariants J1 = I1*I3^-1/3, J2 = I2*I3^-2/3, J3=det F
  std::array<double,3> computeReducedInvariants(const std::array<double,3> &invariants, double deformationGradientDeterminant);

  //! compute the elasticity tensor C = 2*sym(dS/dC). Due to hyperelasticity there are symmetries C_{ijrs} = C_{jirs} and C_{ijrs} = C_{rsij} that leave 21 independent values.
  std::array<double,21> computeElasticityTensor(const std::array<Vec3,3> &rightCauchyGreen,
                                                const std::array<Vec3,3> &inverseRightCauchyGreen,
                                                const std::array<double,3> invariants,
                                                const std::array<double,3> &reducedInvariants);

  //! compute 2nd Piola-Kirchhoff stress tensor S = 2*sym(dPsi/dC)
  std::array<Vec3,3> computePK2Stress(const std::array<Vec3,3> &rightCauchyGreen,
                                      const std::array<Vec3,3> &inverseRightCauchyGreen,
                                      const std::array<double,3> invariants,
                                      const std::array<double,3> &reducedInvariants);

  //! return the index to the elastcity array with 21 distinct entries, for logical component klrs
  int getElasticityEntryNo(int k, int l, int r, int s);

  //! return the entry klrs of the elasticity tensor
  double getElasticityEntry(std::array<double, 21> &elasticity, int k, int l, int r, int s);

  double c1_, c2_;  ///< material constants for mooney-rivlin material
  double kappa_;  ///< bulk modulus

};
*/
};  // namespace

#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_compressible.tpp"
