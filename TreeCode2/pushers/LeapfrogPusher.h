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
#include "../Configuration.h"
#include "../potentials/Potential.h"
#include "../bounds/BoundaryConditions.h"
#include "../Tree.h"

namespace treecode {

namespace pusher {

template <class Vec, class Mat>
class LeapfrogPusher : public Pusher<Vec, Mat> {
public:
	/**
	 * @class LeapfrogPusher "pushers/LeapfrogPusher.h"
	 * @brief A particle pusher using the leapfrog method.
	 */

	/**
	 * @brief Construct a new LeapfrogPusher.
	 * @param conf	Global configuration.
	 * @param bc	Boundary conditions.
	 * @param pot	Potential.
	 */
	LeapfrogPusher(const Configuration<Vec>& conf, const BoundaryConditions<Vec>& bc, const potentials::Potential<Vec,Mat>& pot):
		configuration(conf), boundary(bc), potential(pot){}

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
	void init(std::vector<Particle<Vec>*> parts, Tree<Vec,Mat>& tree, potentials::Precision precision){
		using std::vector;
		//Rebuild the tree, to get the right nodes.
		tree.rebuild();
		//Create interaction list
		vector<Node<Vec,Mat>*> ilist;
		//Loop over all particles, get ilist and then push particle.
		for(Particle<Vec>* p : parts){
			tree.getInteractionList(*p, ilist);
			for(Node<Vec,Mat>* n : ilist){
				push_velocity(*p, *n, configuration.getTimestep()/2, precision);
			}
		}
	}

	/**
	 * @brief Push particle, changing position and velocity.
	 * @param parts		List of particles to push
	 * @param tree		Tree to use to push the particles.
	 * @param precision	Precision to use in force calculation.
	 */
	std::pair<double, double> push_particles(std::vector<Particle<Vec>*> parts, Tree<Vec,Mat>& tree, BoundaryConditions<Vec>& bc,
			potentials::Precision precision){

		using Eigen::Vector2d;

		double ke = 0;
		double pe = 0;


	//	for(Particle* p : parts){
	//		p->updatePosition(p->getVelocity() * 0.001);
	//		std::cout << p->getPosition()[0] << "\t" << p->getPosition()[1] << "\t";
	//	}
	//
	//	for(Particle* p1 : parts){
	//		for(Particle* p2 : parts){
	//			if(p1 == p2)continue;
	//			Vector2d dr = p1->getPosition() - p2->getPosition();
	//			double softened_r = 1.0 / sqrt(dr.squaredNorm() + 0.1*0.1);
	//			Vector2d force = p1->getCharge() * p2->getCharge() * dr * (softened_r*softened_r*softened_r);
	//			p1->updateVelocity(force / p1->getMass() * 0.001);
	//		}
	//	}

		//Push position of each particle

		for(Particle<Vec>* p : parts){
			push_position(*p, configuration.getTimestep());
			bc.particleMoved(p);
		}

		//Now push the velocity
		tree.rebuild();

	//	for(Particle* p : parts){
	//	std::for_each(parts.begin(), parts.end(), [&tree,this,&precision] (Particle* p){
		#pragma omp parallel for reduction(+:ke,pe) schedule(dynamic)
		for(unsigned int i=0;i<parts.size();i++){
			Particle<Vec>* p = parts[i];
			std::vector<Node<Vec,Mat>*> ilist;
			tree.getInteractionList(*p, ilist);
			double initial_vel = p->getVelocity().norm();
			for(Node<Vec,Mat>* n : ilist){
				push_velocity(*p, *n, configuration.getTimestep(), precision);
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
	void push_velocity(Particle<Vec>& particle, const Node<Vec,Mat>& node, double dt, potentials::Precision precision){
		Vec force = potential.getForce(particle, node, precision);
	//	std::cerr << force[0] << "\t" << force[1] << "\t" << force[2] << "\t";
		particle.updateVelocity(force / particle.getMass() * dt);
	}

	/**
	 * @brief Push position of particle, using particle velocity.
	 * @param particle	Particle to push.
	 * @param dt		Timestep.
	 */
	void push_position(Particle<Vec>& particle, double dt){
		particle.updatePosition(particle.getVelocity() * dt);
	}

	const Configuration<Vec>& configuration;
	const BoundaryConditions<Vec>& boundary;
	const potentials::Potential<Vec,Mat>& potential;
};

} /* namespace pushers */
} /* namespace treecode */
#endif /* LEAPFROGPUSHER_H_ */
