///@file

#ifndef MAXWELLDISTRIBUTION_H_
#define MAXWELLDISTRIBUTION_H_

#include <Eigen/Dense>
#include <boost/random/normal_distribution.hpp>
#include <cmath>

#include "Distribution.h"

namespace treecode {
namespace distribution {

/**
 * @class MaxwellDistribution "distributions/MaxwellDistribution.h"
 * @brief Generates vectors according to a Maxwell distribution.
 * @tparam RNG Boost random number generator.
 */
template <class RNG>
class MaxwellDistribution : public VectorDistribution{
public:
	/**
	 * @brief Construct new MaxwellDistribution.
	 * @param c		Global configuration.
	 * @param mass	Mass of particles.
	 * @param temp Temperature of system (as a multiple of the "natural" temperature).
	 */
	MaxwellDistribution(RNG& rng, double mass, double temp, int dims) :
		rng_(rng), mass_(mass), temperature_(temp), D(dims){}

	/**
	 * @brief Get random vectors following a Maxwell distribution.
	 * @param rng	Random number generator.
	 * @return	Vectors following Maxwell Distribution.
	 */
	virtual Eigen::VectorXd getVector() const {
		boost::normal_distribution<double> dist(0.0, sqrt(temperature_ / mass_));
		Eigen::VectorXd v(D);

		for (unsigned int j = 0; j < v.rows(); j++)
			v[j] = dist(rng_);
		return v;
	}

	/**
	 * Unused destructor.
	 */
	virtual ~MaxwellDistribution(){}
private:
	RNG& rng_;
	double mass_, temperature_;
	int D;
};

}
}

#endif /* MAXWELLDISTRIBUTION_H_ */
