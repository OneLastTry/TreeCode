///@file

#ifndef UNIFORMDISTRIBUTION_H_
#define UNIFORMDISTRIBUTION_H_

#include <Eigen/Dense>
#include <boost/random/uniform_01.hpp>

#include "Distribution.h"

#include <iostream>

namespace treecode {
namespace distribution {

/**
 * @class UniformDistribution "distributions/UniformDistribution.h"
 * @brief A class representing a uniform distribution.
 * @tparam RNG Boost random number generator.
 */

template <class RNG>
class UniformDistribution : public VectorDistribution{
public:
	/**
	 * @brief Construct new uniform distribution.
	 * @param min	Minimum points of distribution.
	 * @param max	Maximum points of distribution.
	 */
	UniformDistribution(RNG& rng, const Eigen::VectorXd& min, const Eigen::VectorXd& max) :
		rng_(rng), minimum(min), maximum(max){}

	/**
	 * @brief Generate random vector, with each component somewhere between compoments of min and max.
	 * @param rng	Random number generator.
	 * @return	Randomly distributed vector.
	 */
	virtual Eigen::VectorXd getVector() const  {
		boost::uniform_01<double> dist;
		Eigen::VectorXd v(minimum.rows());
		for (unsigned int j = 0; j < v.rows(); j++) {
			v[j] = dist(rng_) * (maximum[j] - minimum[j]) + minimum[j];
		}
		return v;
	}
	/**
	 * Destructor. Does nothing.
	 */
	virtual ~UniformDistribution(){}
private:
	RNG& rng_;
	Eigen::VectorXd minimum, maximum;
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* UNIFORMDISTRIBUTION_H_ */
