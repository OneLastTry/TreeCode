///@file

#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include <Eigen/Dense>

namespace treecode {
namespace distribution {

/**
 * @class VectorDistribution "distributions/Distribution.h"
 * @brief Base class for distributions returning a vector.
 * @tparam RNG Boost random number generator class (eg, mt19937).
 */
template<class RNG, int D>
class VectorDistribution {
	typedef Eigen::Matrix<double, D, 1> Vec;

public:
	/**
	 * @brief Pure virtual method that must return a vector.
	 * @param rng	Random number generator.
	 * @return Vector in derived class distribution.
	 */
	virtual Vec getVector(RNG& rng) const = 0;
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* DISTRIBUTION_H_ */
