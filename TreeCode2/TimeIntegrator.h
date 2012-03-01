/*
 * TimeIntegrator.h
 *
 *  Created on: 22 Dec 2011
 *      Author: stefans
 */

#ifndef TIMEINTEGRATOR_H_
#define TIMEINTEGRATOR_H_

#include <vector>
#include <string>
#include <fstream>
#include "Particle.h"
#include "Tree.h"
#include "Node.h"
#include "Configuration.h"
#include "bounds/BoundaryConditions.h"
#include "pushers/pusher.h"
#include "potentials/Potential.h"

namespace treecode {

template <class Vec, class Mat>
class TimeIntegrator {
public:
	/**
	 * @class TimeIntegrator
	 * @brief Class for the time integration of particles.
	 *
	 * The purpose of this class is to tie together, with
	 * a supplied particle pusher, various
	 * objects, such as potential, tree and boundary conditions.
	 */

	/**
	 * @brief Instantiate new TimeIntegrator.
	 * @param conf			Configuration.
	 * @param particles		List of particles.
	 * @param tree			Tree.
	 * @param bounds		Boundary conditions.
	 * @param pusher		Particle pusher.
	 */
	TimeIntegrator(const Configuration<Vec>& conf,
			std::vector<Particle<Vec>*>& particles,
			Tree<Vec,Mat>& tree,
			BoundaryConditions<Vec>& bounds,
			pusher::Pusher<Vec,Mat>& pusher):
		bounds_(bounds), particles_(particles), conf_(conf), pusher_(pusher), tree_(tree),
		pos_out_(NULL), vel_out_(NULL), energies_out_(NULL)
	{}

	/**
	 * @brief Kick off the integration.
	 * @param precision		Precision to supply to the potential.
	 * @param output_every	Output every nth step.
	 * @param pos_file		File to output particle positions to (or NULL).
	 * @param vel_file		File to output particle velocities to (or NULL).
	 * @param energy_file	File to output energies to (or NULL).
	 */
	void start(potentials::Precision precision, unsigned int output_every){
		long int num_steps = conf_.getMaxTime() / conf_.getTimestep();
		for(long int i=0;i<num_steps;i++){
			std::pair<double, double> energies = pusher_.push_particles(particles_, tree_, bounds_, precision);
			bounds_.timestepOver();

			std::cout << "Timestep " << i << " of " << num_steps << " complete (" << ((float)i/(float)num_steps)*100 << "%)" << std::endl;
			if( (i%output_every) == 0){
				if(energies_out_ != NULL){
					(*energies_out_) << energies.first << "\t" << energies.second << std::endl;
				}
				if(pos_out_ != NULL || vel_out_ != NULL){
					BOOST_FOREACH(Particle<Vec>* p, particles_){
						if(pos_out_ != NULL){
							for(int i=0;i<p->getPosition().rows();i++)
								(*pos_out_) << p->getPosition()[i] << "\t";
						}if(vel_out_ != NULL){
							for(int i=0;i<p->getVelocity().rows();i++)
								(*vel_out_) << p->getVelocity()[i] << "\t";
						}
					}
				}
			}
		}
	}

	/**
	 * @brief Destructor.
	 */
	~TimeIntegrator() {
		delete pos_out_;
		delete vel_out_;
		delete energies_out_;
	}

	/**
	 * @brief Set file to output particle positions to
	 * @param filename File name
	 */
	void setPositionOutputFile(const char* filename){
		pos_out_ = new std::ofstream(filename);
	}

	/**
	 * @brief Set file to output particle velocities to
	 * @param filename File name
	 */
	void setVelocityOutputFile(const char* filename){
		vel_out_ = new std::ofstream(filename);
	}

	/**
	 * @brief Set file to output system energies to
	 * @param filename File name
	 */
	void setEnergyOutputFile(const char* filename){
		energies_out_ = new std::ofstream(filename);
	}

private:
	BoundaryConditions<Vec>& bounds_;
	std::vector<Particle<Vec>*>& particles_;
	const Configuration<Vec>& conf_;
	pusher::Pusher<Vec,Mat>& pusher_;
	Tree<Vec,Mat>& tree_;

	std::ofstream *pos_out_, *vel_out_, *energies_out_;
};

} /* namespace treecode */
#endif /* TIMEINTEGRATOR_H_ */
