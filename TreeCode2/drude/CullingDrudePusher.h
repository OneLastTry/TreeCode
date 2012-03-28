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

template<int D, class RNG>
class CullingDrudePusher : public DrudePusher<D>{
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;

public:
	CullingDrudePusher(const Configuration<Vec>& config,
			BoundaryConditions<D>& bounds,
			const distribution::VectorDistribution<RNG, D>& vel_dist,
			DrudeMAC<D>& mac,
			RNG& rng):
				DrudePusher<D>(config, bounds, mac),
				vel_dist_(vel_dist), rng_(rng){}

	std::pair<double, double> push_particles(std::vector<Particle<D>*> parts, Tree<D>& tree, BoundaryConditions<D>& bc,
			potentials::Precision prec, const AcceptanceCriterion<D>& mac){

		//If next particle step would take it across the entire system, reset velocity vector
		typedef Particle<D> part_t;
		BOOST_FOREACH(part_t* p, parts){
			if(p->getVelocity().norm() * DrudePusher<D>::config_.getTimestep() > DrudePusher<D>::bounds_.getSize()){
				p->setVelocity(vel_dist_.getVector(rng_));
				std::cout << "Culling particle" << std::endl;
			}
		}

		return DrudePusher<D>::push_particles(parts, tree, bc, prec, mac);
	}
protected:
	const distribution::VectorDistribution<RNG, D>& vel_dist_;
	RNG& rng_;
};

}
}

#endif /* CULLINGDRUDEPUSHER_H_ */
