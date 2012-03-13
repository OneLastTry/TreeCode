/*
 * DrudePusher.h
 *
 *  Created on: 10 Mar 2012
 *      Author: stefans
 */

#ifndef DRUDEPUSHER_H_
#define DRUDEPUSHER_H_

#include "../pushers/pusher.h"
#include "../Configuration.h"
#include "../Node.h"
#include "../Particle.h"
#include "BigParticle.h"
#include <vector>
#include <boost/foreach.hpp>
#include <cmath>

namespace treecode{
namespace pusher{

template<class Vec, class Mat>
class DrudePusher : public Pusher<Vec, Mat>{
public:
	DrudePusher(const Configuration<Vec>& config, BoundaryConditions<Vec>& bounds):config_(config), bounds_(bounds){}

	std::pair<double, double> push_particles(std::vector<Particle<Vec>*> parts, Tree<Vec,Mat>& tree,
			BoundaryConditions<Vec>& bc,
			potentials::Precision prec, const AcceptanceCriterion<Vec, Mat>& mac){

		std::pair<double, double> energies(0,0);
		BOOST_FOREACH(Particle<Vec>* p, parts){
				std::pair<double, double> p_energy = push_particle_recursive(p, tree, bc, prec, mac, config_.getTimestep());
				energies.first += p_energy.first;
				energies.second += p_energy.second;
		}
		return energies;
	}

private:
	std::pair<double, double> push_particle_recursive(Particle<Vec>* p, Tree<Vec,Mat>& tree,
			BoundaryConditions<Vec>& bc,
				potentials::Precision prec, const AcceptanceCriterion<Vec, Mat>& mac, double dt){
		typedef std::vector<Node<Vec,Mat>* > interaction_list;

		tree.rebuild();

		double ke = 0;
		interaction_list ilist;
		tree.getInteractionList(*p, ilist, mac);

		if(ilist.size() > 0){
			double distance_to_travel = p->getVelocity().norm() * config_.getTimestep();
			for(typename interaction_list::iterator it = ilist.begin(); it < ilist.end(); it++){
				Node<Vec,Mat>* n = *it;
				BigParticle<Vec>* interacting_particle = dynamic_cast<BigParticle<Vec>*>(n->getParticles().front());
				double distance_to_edge = distanceToIntersection(p, interacting_particle);
				//Move to edge of particle
				p->updatePosition(p->getVelocity() * dt * (distance_to_edge/distance_to_travel));
				//Reflect the velocity vector
				Vec disp_vec = p->getPosition() - interacting_particle->getPosition();
				disp_vec /= disp_vec.norm();
				Mat reflection_matrix = Mat::Identity() - 2.0 * disp_vec * disp_vec.transpose();
				p->setVelocity(reflection_matrix * p->getVelocity());
				//Now move forward again
				push_particle_recursive(p, tree, bc, prec, mac, dt * (1 - distance_to_edge / distance_to_travel));
			}
		}else{
			p->updatePosition(p->getVelocity() * dt);
		}
		ke += 0.5 * p->getMass() * p->getVelocity().squaredNorm();
		bounds_.particleMoved(p);

		return std::pair<double, double>(ke, 0);
	}

	double distanceToIntersection(Particle<Vec>* p, BigParticle<Vec>* bp){
		Vec l = p->getVelocity() / p->getVelocity().norm();
		Vec c = bp->getPosition() - p->getPosition();
		return l.dot(c) - sqrt(l.dot(c) * l.dot(c) - c.squaredNorm() + bp->getRadius()*bp->getRadius());
	}

protected:
	const Configuration<Vec>& config_;
	BoundaryConditions<Vec>& bounds_;
};

}
}

#endif /* DRUDEPUSHER_H_ */
