///@file

#ifndef SINUSOIDALDISTRIBUTION_H_
#define SINUSOIDALDISTRIBUTION_H_

#include "Distribution.h"
#include "UniformDistribution.h"
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
template <class RNG>
class SinusoidalDistribution : public VectorDistribution{

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
			RNG& rng,
			unsigned int dimension,
			const Eigen::VectorXd& min,
			const Eigen::VectorXd& max,
			double wavelengths, double phase):
				rng_(rng),
				dimension_(dimension),
				min_(min), max_(max),
				wavelengths_(wavelengths), phase_(phase),
				real_dist_(0, 2){
		uniform_dist_ = new UniformDistribution<RNG>(rng, min, max);
	}

	virtual Eigen::VectorXd getVector() const{
		Eigen::VectorXd v(min_.rows());
		//Wavelength is the size of the system / the number of wavelengths
		double wavelength = (max_[dimension_] - min_[dimension_]) / wavelengths_;
		double k = 2.0 * M_PI / wavelength;	//Wavenumber
		double u;
		do{
			u = real_dist_(rng_);
			v = uniform_dist_->getVector();
		}while(sin(v[dimension_]*k + phase_)+1 < u);
		return v;
	}

	virtual ~SinusoidalDistribution(){
		delete uniform_dist_;
	}
private:
	RNG& rng_;
	UniformDistribution<RNG>* uniform_dist_;
	unsigned int dimension_;
	Eigen::VectorXd min_, max_;
	double wavelengths_, phase_;
	boost::random::uniform_real_distribution<double> real_dist_;
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* SINUSOIDALDISTRIBUTION_H_ */
