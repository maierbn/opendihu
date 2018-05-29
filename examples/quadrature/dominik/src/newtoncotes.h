#pragma once

//#include <Python.h>  // has to be the first included header
#include <array>
#include <type_traits>

#include "integrator.h"

namespace Quadrature
{

	/** Implements Gauss quadrature with NumberGaussPoints. It is capable of exactly integrating polynomials of degree 2*NumberGaussPoints-1.
	*/
	template<unsigned int NumberNCPoints>
	class NC : public Quadrature
	{
	public:
		typedef NC<NumberNCPoints> HighOrderQuadrature;   ///< this defines the own class, to be able to generalize code to mixed quadrature

																//! return the number of evaluations that are needed for a 1D quadrature
		static constexpr int numberEvaluationsD();
		static constexpr int numberEvaluations = numNodes;

		//! return the sampling points, i.e. newton-cotes points that are needed for the quadrature. The list may not be in ascending order, but the order matches the order required in integrate
		static std::array<double, NumberNCPoints> samplingPoints();

		//! Compute the integral from evaluations at the gauss points.
		//! If a std::array is given for ValueType, compute separate integrals for each component with the same gauss points for all.
		template<typename ValueType>
		//static ValueType computeIntegral(const typename std::array<ValueType,numberEvaluations()>::const_iterator evaluations);
		static ValueType computeIntegral(const typename std::array<ValueType, numberEvaluations>::const_iterator evaluations);

		//! Compute the integral from evaluations at the gauss points.
		//! If a std::array is given for ValueType, compute separate integrals for each component with the same gauss points for all.
		template<typename ValueType>
		//tatic ValueType computeIntegral(const typename std::array<ValueType,numberEvaluations()> &evaluations);
		static ValueType computeIntegral(const typename std::array<ValueType, numberEvaluations> &evaluations);
	};

};

#include "newtoncotes.tpp"