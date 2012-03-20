///@file

#ifndef UNIFORMDISTRIBUTION_H_
#define UNIFORMDISTRIBUTION_H_

#include <Eigen/Dense>
#include <boost/random/uniform_01.hpp>

#include "Distribution.h"

namespace treecode {
namespace distribution {

/**
 * @class UniformDistribution "distributions/UniformDistribution.h"
 * @brief A class representing a uniform distribution.
 * @tparam RNG Boost random number generator.
 */

template <class RNG, class Vec>
class UniformDistribution : public VectorDistribution<RNG, Vec>{
public:
	/**
	 * @brief Construct new uniform distribution.
	 * @param min	Minimum points of distribution.
	 * @param max	Maximum points of distribution.
	 */
	UniformDistribution(const Vec& min, const Vec& max) : minimum(min), maximum(max){}

	/**
	 * @brief Generate random vector, with each component somewhere between compoments of min and max.
	 * @param rng	Random number generator.
	 * @return	Randomly distributed vector.
	 */
	virtual Vec getVector(RNG& rng) const  {
		boost::uniform_01<double> dist;
		Vec v;
		for (unsigned int j = 0; j < v.rows(); j++) {
			v[j] = dist(rng) * (maximum[j] - minimum[j]) + minimum[j];
		}
		return v;
	}
	/**
	 * Destructor. Does nothing.
	 */
	virtual ~UniformDistribution(){}
private:
	Vec minimum, maximum;
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* UNIFORMDISTRIBUTION_H_ */
