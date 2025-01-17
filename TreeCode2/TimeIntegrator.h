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
#include "bounds/BoundaryConditions.h"
#include "pushers/pusher.h"
#include "potentials/Potential.h"
#include "io/ParticleTracker.h"

namespace treecode {

template <int D>
class TimeIntegrator {
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;

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
	 * @param particles		List of particles.
	 * @param tree			Tree.
	 * @param bounds		Boundary conditions.
	 * @param pusher		Particle pusher.
	 */
	TimeIntegrator(
			double timestep, double max_time,
			std::vector<Particle<D>*>& particles,
			Tree<D>& tree,
			BoundaryConditions<D>& bounds,
			pusher::Pusher<D>& pusher,
			const AcceptanceCriterion<D>& mac):
		dt_(timestep), max_time_(max_time),
		bounds_(bounds), particles_(particles), pusher_(pusher),
		tree_(tree), energies_out_(NULL), mac_(mac)
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
		long int num_steps = max_time_/dt_;
		for(long int i=0;i<num_steps;i++){

			std::pair<double, double> energies = pusher_.push_particles(particles_, tree_, bounds_, precision, mac_);
			bounds_.timestepOver();

			std::cout << "\rTimestep " << i << " of " << num_steps << " complete (" << ((float)i/(float)num_steps)*100 << "%)";
			if( (i%output_every) == 0){
				if(energies_out_ != NULL)
					(*energies_out_) << (i * dt_) << "\t" << energies.first << "\t" << energies.second << std::endl;

				typedef output::ParticleTracker<D> track;
				BOOST_FOREACH(track* t, particle_trackers_)
					t->output();
			}
		}
		std::cout << std::endl;
	}

	/**
	 * @brief Destructor.
	 */
	~TimeIntegrator() {
		delete energies_out_;
	}

	/**
	 * @brief Set file to output system energies to
	 * @param filename File name
	 */
	void setEnergyOutputFile(const char* filename){
		energies_out_ = new std::ofstream(filename);
	}

	/**
	 * @brief Add particle tracker.
	 * @param tracker Tracker.
	 */
	void addParticleTracker(output::ParticleTracker<D>* tracker){
		particle_trackers_.push_back(tracker);
	}

	template<class InputIterator>
	void addParticleTrackers(InputIterator begin, InputIterator end){
		particle_trackers_.insert(particle_trackers_.end(), begin, end);
	}

private:
	double dt_, max_time_;
	BoundaryConditions<D>& bounds_;
	std::vector<Particle<D>*>& particles_;
	pusher::Pusher<D>& pusher_;
	Tree<D>& tree_;

	std::ofstream *energies_out_;
	const AcceptanceCriterion<D>& mac_;

	std::vector<output::ParticleTracker<D>* > particle_trackers_;
};

} /* namespace treecode */
#endif /* TIMEINTEGRATOR_H_ */
