///@file

#ifndef CONSTDISTRIBUTION_H_
#define CONSTDISTRIBUTION_H_

#include "Distribution.h"
#include "../Configuration.h"

#include <Eigen/Dense>

namespace treecode {
namespace distribution {

/**
 * @class ConstDistribution "distributions/ConstDistribution.h"
 * @brief Distribution class that returns a constant vector.
 * @tparam Boost random number generator.
 */
template <class RNG, class Vec>
class ConstDistribution : public VectorDistribution<RNG, Vec>{
public:
	/**
	 * @brief Construct new ConstDistribution.
	 * @param _constant	Vector that will always be returned by ConstDistribution::getVector().
	 */
	ConstDistribution(const Vec& _constant) : constant(_constant) {}

	/**
	 * @brief Return constant vector, as initialised by constructor.
	 * @param rng	(Unnecessary) random number generator.
	 * @return Constant vector.
	 */
	virtual Vec getVector(RNG& rng) const {
		return constant;
	}
	/**
	 * Destructor. Unused.
	 */
	virtual ~ConstDistribution(){}
private:
	const Vec constant;
};

}
}

#endif /* CONSTDISTRIBUTION_H_ */
