#include "newtoncotes.h"

#include <array>
#include <cmath>

namespace Quadrature
{
	// 1 NC point --> rectangle
	template<>
	std::array<double, 1> NC<1>::
		samplingPoints()
	{
		return std::array<double, 1>{
			0.5
		};
	}

	// 2 NC points --> trapezoidal
	template<>
	std::array<double, 2> NC<2>::
		samplingPoints()
	{
		return std::array<double, 2>{
			0.,
			1.
		};
	}

	// 3 NC points --> simpson
	template<>
	std::array<double, 3> NC<3>::
		samplingPoints()
	{
		return std::array<double, 3>{
			0., 
				0.5, 
				1.
		};
	}

	// 4 NC points --> 3/8
	template<>
	std::array<double, 4> NC<4>::
		samplingPoints()
	{
		return std::array<double, 4>{
			0., 
				(1. / 3.), 
				(2. / 3.), 
				1.
		};
	}

	// 5 NC points
	template<>
	std::array<double, 5> NC<5>::
		samplingPoints()
	{
		return std::array<double, 5>{
			0.,
				1./4.,
				2./4.,
				3./4.,
				1.
		};
	}

	// 6 NC points
	template<>
	std::array<double, 6> NC<6>::
		samplingPoints()
	{
		return std::array<double, 6>{
			0.,
				1./5.,
				2./5.,
				3./5.,
				4./5.,
				1.,
		};
	}

	// 7 NC points
	template<>
	std::array<double, 7> NC<7>::
		samplingPoints()
	{
		return std::array<double, 7>{
			0.,
				1./6.,
				2./6.,
				3./6.,
				4./6.,
				5./6.,
				1.,
		};
	}

	// 8 NC points
	template<>
	std::array<double, 8> NC<8>::
		samplingPoints()
	{
		return std::array<double, 8>{
			0.,
				1. / 7.,
				2. / 7.,
				3. / 7.,
				4. / 7.,
				5. / 7.,
				6. / 7.,
				1.
		};
	}

};