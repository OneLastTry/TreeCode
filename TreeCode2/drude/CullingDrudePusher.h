/*
 * DrudePusher.h
 *
 *  Created on: 10 Mar 2012
 *      Author: stefans
 */

#ifndef CULLINGDRUDEPUSHER_H_
#define CULLINGDRUDEPUSHER_H_

#include "DrudePusher.h"
#include "../distributions/Distribution.h"
#include "../Particle.h"
#include <boost/foreach.hpp>

#include <iostream>

namespace treecode{
namespace pusher{

template<class Vec, class Mat, class RNG>
class CullingDrudePusher : public DrudePusher<Vec, Mat>{
public:
	CullingDrudePusher(const Configuration<Vec>& config,
			BoundaryConditions<Vec,Mat>& bounds,
			const distribution::VectorDistribution<RNG, Vec>& vel_dist,
			DrudeMAC<Vec,Mat>& mac,
			RNG& rng):
				DrudePusher<Vec, Mat>(config, bounds, mac),
				vel_dist_(vel_dist), rng_(rng){}

	std::pair<double, double> push_particles(std::vector<Particle<Vec,Mat>*> parts, Tree<Vec,Mat>& tree, BoundaryConditions<Vec,Mat>& bc,
			potentials::Precision prec, const AcceptanceCriterion<Vec, Mat>& mac){

		//If next particle step would take it across the entire system, reset velocity vector
		typedef Particle<Vec,Mat> part_t;
		BOOST_FOREACH(part_t* p, parts){
			if(p->getVelocity().norm() * DrudePusher<Vec,Mat>::config_.getTimestep() > DrudePusher<Vec,Mat>::bounds_.getSize()){
				p->setVelocity(vel_dist_.getVector(rng_));
				std::cout << "Culling particle" << std::endl;
			}
		}

		return DrudePusher<Vec, Mat>::push_particles(parts, tree, bc, prec, mac);
	}
protected:
	const distribution::VectorDistribution<RNG, Vec>& vel_dist_;
	RNG& rng_;
};

}
}

#endif /* CULLINGDRUDEPUSHER_H_ */
