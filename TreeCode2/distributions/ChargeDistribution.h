///@file

#ifndef CHARGEDISTRIBUTION_H_
#define CHARGEDISTRIBUTION_H_

namespace treecode {
namespace distribution {

/**
 * @class ChargeDistribution "distributions/ChargeDistribution.h"
 * @brief Base class for distributions of particle charge.
 * @tparam RNG Boost random number generator.
 */
class ChargeDistribution {
public:
	/**
	 * @brief Pure virtual function that should return charge.
	 * @param rng	Random number generator.
	 * @return Charge corresponding to derived class distribution.
	 */
	virtual int getCharge() const = 0;
	virtual ~ChargeDistribution(){}
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* CHARGEDISTRIBUTION_H_ */
