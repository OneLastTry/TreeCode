///@file

#ifndef MAXWELLDISTRIBUTION_H_
#define MAXWELLDISTRIBUTION_H_

#include <Eigen/Dense>
#include <boost/random/normal_distribution.hpp>
#include <cmath>

#include "Distribution.h"
#include "../Configuration.h"

namespace treecode {
namespace distribution {

/**
 * @class MaxwellDistribution "distributions/MaxwellDistribution.h"
 * @brief Generates vectors according to a Maxwell distribution.
 * @tparam RNG Boost random number generator.
 */
template <class RNG, class Vec>
class MaxwellDistribution : public VectorDistribution<RNG, Vec>{
public:
	/**
	 * @brief Construct new MaxwellDistribution.
	 * @param c		Global configuration.
	 * @param mass	Mass of particles.
	 * @param temp Temperature of system (as a multiple of the "natural" temperature).
	 */
	MaxwellDistribution(double mass, double temp) : mass_(mass), temperature_(temp){}

	/**
	 * @brief Get random vectors following a Maxwell distribution.
	 * @param rng	Random number generator.
	 * @return	Vectors following Maxwell Distribution.
	 */
	virtual Vec getVector(RNG& rng) const {
		boost::normal_distribution<double> dist(0.0, sqrt(temperature_ / mass_));
		Vec v;

		for (unsigned int j = 0; j < v.rows(); j++)
			v[j] = dist(rng);
		return v;
	}

	/**
	 * Unused destructor.
	 */
	virtual ~MaxwellDistribution(){}
private:
	double mass_, temperature_;
};

}
}

#endif /* MAXWELLDISTRIBUTION_H_ */
