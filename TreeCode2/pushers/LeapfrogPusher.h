/*
 * LeapfrogPusher.h
 *
 *  Created on: 7 Dec 2011
 *      Author: stefans
 */

#ifndef LEAPFROGPUSHER_H_
#define LEAPFROGPUSHER_H_

#include "pusher.h"
#include "../potentials/Potential.h"

#include "../Particle.h"
#include "../Node.h"
#include "../potentials/Potential.h"
#include "../bounds/BoundaryConditions.h"
#include "../Tree.h"
#include "../macs/AcceptanceCriterion.h"

namespace treecode {

namespace pusher {

template <int D>
class LeapfrogPusher : public Pusher<D> {
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;
public:
	/**
	 * @class LeapfrogPusher "pushers/LeapfrogPusher.h"
	 * @brief A particle pusher using the leapfrog method.
	 */

	/**
	 * @brief Construct a new LeapfrogPusher.
	 * @param bc	Boundary conditions.
	 * @param pot	Potential.
	 */
	LeapfrogPusher(double timestep, const BoundaryConditions<D>& bc, const potentials::Potential<D>& pot):
		dt_(timestep), boundary(bc), potential(pot){}

	/**
	 * @brief Initialise leapfrog pusher.
	 *
	 * This causes the pusher to push the velocity of all particles
	 * by half a timestep (using a simple Euler step) in order to
	 * get the required half step lag.
	 *
	 * @param parts		List of particles that must be initialised.
	 * @param tree		Tree object, containg all particles.
	 * @param precision	Precision to use in force calculation.
	 */
	void init(std::vector<Particle<D>*> parts, Tree<D>& tree, potentials::Precision precision,
			const AcceptanceCriterion<D>& mac){
		using std::vector;
		typedef Node<D> Node;
		typedef vector<Node*> interaction_list;
		typedef Particle<D> part_t;

		//Rebuild the tree, to get the right nodes.
		tree.rebuild();
		//Loop over all particles, get ilist and then push particle.
#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<parts.size();i++){
			interaction_list ilist;
			part_t* p = parts[i];
			tree.getInteractionList(*p, ilist, mac);
			for(typename interaction_list::iterator it = ilist.begin(); it < ilist.end(); it++){
				Node* n = *it;
				push_velocity(*p, *n, dt_/2, precision);
			}
		}
	}

	/**
	 * @brief Push particle, changing position and velocity.
	 * @param parts		List of particles to push
	 * @param tree		Tree to use to push the particles.
	 * @param precision	Precision to use in force calculation.
	 */
	std::pair<double, double> push_particles(
			std::vector<Particle<D>*> parts, Tree<D>& tree, BoundaryConditions<D>& bc,
			potentials::Precision precision,
			const AcceptanceCriterion<D>& mac){

		typedef Node<D> Node;
		typedef std::vector<Node*> interaction_list;
		typedef Particle<D> part_t;

		using Eigen::Vector2d;

		double ke = 0;
		double pe = 0;

		//Push position of each particle

		BOOST_FOREACH(part_t* p, parts){
			push_position(*p, dt_);
			bc.particleMoved(p);
		}

		//Now push the velocity
		tree.rebuild();

		#pragma omp parallel for reduction(+:ke,pe) schedule(dynamic)
		for(unsigned int i=0;i<parts.size();i++){
			Particle<D>* p = parts[i];
			interaction_list ilist;
			tree.getInteractionList(*p, ilist, mac);
			double initial_vel = p->getVelocity().norm();
			for(typename interaction_list::iterator it = ilist.begin(); it < ilist.end(); it++){
				Node* n = *it;
				push_velocity(*p, *n, dt_, precision);
				pe += 0.5 * p->getCharge() * potential.getPotential(*p, *n, precision);
			}
			double new_vel = p->getVelocity().norm();
			double mean_vel = (initial_vel + new_vel)/2;
			ke += 0.5 * p->getMass() * mean_vel*mean_vel;
		}
		return std::pair<double, double>(ke, pe);
	}

	/**
	 * Destructor. Does nothing.
	 */
	~LeapfrogPusher() {}

private:
	/**
	 * @brief Push particle velocity, getting force from Potential object.
	 * @param particle	Particle to push.
	 * @param node		Node the particle is interacting with.
	 * @param dt		Timestep.
	 * @param precision	Precision to use in force calculation.
	 */
	void push_velocity(Particle<D>& particle, const Node<D>& node, double dt, potentials::Precision precision){
		Vec force = potential.getForce(particle, node, precision);
	//	std::cerr << force[0] << "\t" << force[1] << "\t" << force[2] << "\t";
		particle.updateVelocity(force / particle.getMass() * dt);
	}

	/**
	 * @brief Push position of particle, using particle velocity.
	 * @param particle	Particle to push.
	 * @param dt		Timestep.
	 */
	void push_position(Particle<D>& particle, double dt){
		particle.updatePosition(particle.getVelocity() * dt);
	}

	double dt_;
	const BoundaryConditions<D>& boundary;
	const potentials::Potential<D>& potential;
};

} /* namespace pushers */
} /* namespace treecode */
#endif /* LEAPFROGPUSHER_H_ */
