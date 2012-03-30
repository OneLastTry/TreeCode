///@file

#ifndef CYLINDRICALDISTRIBUTION_H_
#define CYLINDRICALDISTRIBUTION_H_

#include "Distribution.h"
#include "SphericalDistribution.h"
#include <Eigen/Dense>
#include <boost/random/uniform_real_distribution.hpp>

namespace treecode {
namespace distribution {

/**
 * @class CylindricalDistribution "distributions/CylindricalDistribution.h"
 * @brief Generates a 3D cylinder of random points.
 * @tparam RNG Random number generator.
 */
template <class RNG>
class CylindricalDistribution : public VectorDistribution{
public:
	/**
	 * @brief Instantiate a new cylindrical distribution with the specified parameters.
	 * @param config		Global configuration.
	 * @param bottom_centre	Centre of the bottom (circular) face of the cylinder.
	 * @param radius		Radius of the cylinder ends.
	 * @param height		Height of the cylinder.
	 */
	CylindricalDistribution(
			RNG& rng,
			const Eigen::Vector2d& bottom_centre,
			double radius,
			double height
	):
	bottom_centre_(bottom_centre),
	radius_(radius), height_(height),
	height_dist_(0, height){
		//Create a circular distribution for the radial part.
		circ_dist_ = new SphericalDistribution<RNG>(rng, 2, bottom_centre, radius);
	}

	/**
	 * @brief Find a random point inside a cylinder.
	 * @param rng	Random number generator.
	 * @return		Random vector in cylinder.
	 */
	Eigen::VectorXd getVector() const{
		Eigen::Vector3d v;
		Eigen::Vector2d circ_vec = circ_dist_->getVector();
		v[0] = circ_vec[0];
		v[1] = circ_vec[1];
		v[2] = height_dist_(rng_);
		return v;
	}

	///Destructor
	virtual ~CylindricalDistribution(){
		delete circ_dist_;
	}
private:
	RNG& rng_;
	SphericalDistribution<RNG>* circ_dist_;
	const Eigen::Vector2d bottom_centre_;
	double radius_, height_;
	boost::random::uniform_real_distribution<double> height_dist_;
};

} /* namespace distribution */
} /* namespace treecode */
#endif /* CYLINDRICALDISTRIBUTION_H_ */
