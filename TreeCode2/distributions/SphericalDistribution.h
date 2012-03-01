///@file

#ifndef SPHERICALDISTRIBUTION_H_
#define SPHERICALDISTRIBUTION_H_

#include "Distribution.h"
#include <boost/random/uniform_real_distribution.hpp>
#include <Eigen/Dense>

#include <iostream>

namespace treecode {

namespace distribution {

/**
 * @class SphericalDistribution "distributions/SphericalDistribution.h"
 * @brief Class used for generating points within an n-sphere.
 *
 * This should be used to generate points in a circle, or sphere. Note that
 * it uses a simple rejector method, and so is not very efficient. I see no
 * problem with that, as it is only used at the start of a simulation, and
 * generalises nicely to n dimensions.
 *
 * @tparam RNG Random number generator.
 */
template <class RNG, class Vec>
class SphericalDistribution : public VectorDistribution<RNG, Vec>{
public:
	/**
	 * @brief Constructor for a distribution that supplies points within an n-sphere.
	 * @param dims		Dimensions.
	 * @param centre	Centre point of n-sphere.
	 * @param radius	Radius of n-sphere.
	 */
	SphericalDistribution(unsigned int dims, const Vec& centre, double radius):
		dims_(dims), centre_(centre), radius_(radius), dist_(-radius, radius){}

	/**
	 * @brief Generate points within an n-dimensional sphere.
	 *
	 * This is a simple method, that generates points within an n-cube,
	 * and rejects any that have radius less than the specified radius.
	 *
	 * @param rng	Random number generator.
	 * @return		Random point within the specified n-sphere.
	 */
	Vec getVector(RNG& rng) const{
		using boost::random::uniform_real_distribution;
		Vec v(dims_);
		do{
			for(unsigned int i = 0; i < dims_; i++)
				v[i] = dist_(rng);
		}while(v.norm() > radius_);	//Reject any in corners of cube/square/hypercube
		return v + centre_;	//Shift up by centre
	}

	///Destrutor that does nothing.
	virtual ~SphericalDistribution(){}
private:
	unsigned int dims_;
	Vec centre_;
	double radius_;
	boost::random::uniform_real_distribution<double> dist_;
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* SPHERICALDISTRIBUTION_H_ */
