///@file

#ifndef UNIFORMCHARGEDISTRIBUTION_H_
#define UNIFORMCHARGEDISTRIBUTION_H_

#include "ChargeDistribution.h"

namespace treecode {
namespace distribution {

/**
 * @class ConstantChargeDistribution "distributions/ConstantChargeDistribution.h"
 * @brief A very simple class that returns a constant charge.
 * @tparam RNG Boost random number generator.
 */

template <class RNG>
class ConstantChargeDistribution : public ChargeDistribution<RNG>{
public:
	/**
	 * @brief Generates a new constant charge distribution.
	 * @param _charge Charge to always return.
	 */
	ConstantChargeDistribution(int _charge) : charge(_charge){}

	/**
	 * @brief Returns a constant charge.
	 * @param rng Unnecessary random number generator.
	 * @return Charge.
	 */
	virtual int getCharge(RNG& rng) const  {return charge;}
	///Destructor, does nothing.
	virtual ~ConstantChargeDistribution(){}
private:
	int charge;
};

} /* namespace distribution */
} /* namespace treecode */

#endif /* UNIFORMCHARGEDISTRIBUTION_H_ */
