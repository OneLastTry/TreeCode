///@file

#ifndef CONSTDISTRIBUTION_H_
#define CONSTDISTRIBUTION_H_

#include "Distribution.h"

#include <Eigen/Dense>

namespace treecode {
namespace distribution {

/**
 * @class ConstDistribution "distributions/ConstDistribution.h"
 * @brief Distribution class that returns a constant vector.
 * @tparam Boost random number generator.
 */
class ConstDistribution : public VectorDistribution{

public:
	/**
	 * @brief Construct new ConstDistribution.
	 * @param _constant	Vector that will always be returned by ConstDistribution::getVector().
	 */
	ConstDistribution(const Eigen::VectorXd& _constant) : constant(_constant) {}

	/**
	 * @brief Return constant vector, as initialised by constructor.
	 * @param rng	(Unnecessary) random number generator.
	 * @return Constant vector.
	 */
	virtual Eigen::VectorXd getVector() const {
		return constant;
	}
	/**
	 * Destructor. Unused.
	 */
	virtual ~ConstDistribution(){}
private:
	const Eigen::VectorXd constant;
};

}
}

#endif /* CONSTDISTRIBUTION_H_ */
