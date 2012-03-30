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
class VectorDistribution {
public:
	/**
	 * @brief Pure virtual method that must return a vector.
	 * @param rng	Random number generator.
	 * @return Vector in derived class distribution.
	 */
	virtual Eigen::VectorXd getVector() const = 0;
	virtual ~VectorDistribution(){}
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* DISTRIBUTION_H_ */
