///@file

#ifndef SINUSOIDALDISTRIBUTION_H_
#define SINUSOIDALDISTRIBUTION_H_

#include "Distribution.h"
#include "UniformDistribution.h"
#include "../Configuration.h"
#include <cmath>
#include <boost/random/uniform_real_distribution.hpp>
#include <Eigen/Dense>

#include <iostream>

namespace treecode {
namespace distribution {


/**
 * @class SinusoidalDistribution "distributions/SinusoidalDistribution.h"
 * @brief Generates particles with sinusoidal density in one dimension, uniform in others.
 * @tparam RNG Random number generator class
 */
template <class RNG, class Vec>
class SinusoidalDistribution : public VectorDistribution<RNG, Vec>{
public:
	/**
	 * @brief Construct distribution objects that returns vectors sinusoidally distributed in a given dimension.
	 * @param config	Global configuration.
	 * @param dimension	The direction in which to sinusoidally distribute particles (uniform in other directions).
	 * @param min		Minimum position to place particles.
	 * @param max		Maximum position to place particles.
	 * @param wavelengths Number of wavelengths to include.
	 * @param phase		Phase of generator.
	 */
	SinusoidalDistribution(
			unsigned int dimension,
			const Vec min,
			const Vec max,
			double wavelengths, double phase):
				dimension_(dimension),
				min_(min), max_(max),
				wavelengths_(wavelengths), phase_(phase),
				real_dist_(0, 2){
		uniform_dist_ = new UniformDistribution<RNG, Vec>(min, max);
	}

	virtual Vec getVector(RNG& rng) const{
		Vec v;
		//Wavelength is the size of the system / the number of wavelengths
		double wavelength = (max_[dimension_] - min_[dimension_]) / wavelengths_;
		double k = 2.0 * M_PI / wavelength;	//Wavenumber
		double u;
		do{
			u = real_dist_(rng);
			v = uniform_dist_->getVector(rng);
		}while(sin(v[dimension_]*k)+1 < u);
		return v;
	}

	virtual ~SinusoidalDistribution(){
		delete uniform_dist_;
	}
private:
	UniformDistribution<RNG, Vec>* uniform_dist_;
	unsigned int dimension_;
	Vec min_, max_;
	double wavelengths_, phase_;
	boost::random::uniform_real_distribution<double> real_dist_;
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* SINUSOIDALDISTRIBUTION_H_ */
