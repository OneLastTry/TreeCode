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
template <class RNG>
class SphericalDistribution : public VectorDistribution{

public:
	/**
	 * @brief Constructor for a distribution that supplies points within an n-sphere.
	 * @param rng	Random number generator.
	 * @param dims		Dimensions.
	 * @param centre	Centre point of n-sphere.
	 * @param radius	Radius of n-sphere.
	 */
	SphericalDistribution(
			RNG& rng,
			unsigned int dims,
			const Eigen::VectorXd& centre,
			double radius):
				rng_(rng), dims_(dims), centre_(centre), radius_(radius), dist_(-radius, radius){}

	/**
	 * @brief Generate points within an n-dimensional sphere.
	 *
	 * This is a simple method, that generates points within an n-cube,
	 * and rejects any that have radius less than the specified radius.
	 *
	 * @return		Random point within the specified n-sphere.
	 */
	Eigen::VectorXd getVector() const{
		Eigen::VectorXd v(dims_);
		do{
			for(unsigned int i = 0; i < dims_; i++)
				v[i] = dist_(rng_);
		}while(v.norm() > radius_);	//Reject any in corners of cube/square/hypercube
		return v + centre_;	//Shift up by centre
	}

	///Destrutor that does nothing.
	virtual ~SphericalDistribution(){}
private:
	RNG& rng_;
	unsigned int dims_;
	Eigen::VectorXd centre_;
	double radius_;
	boost::random::uniform_real_distribution<double> dist_;
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* SPHERICALDISTRIBUTION_H_ */
