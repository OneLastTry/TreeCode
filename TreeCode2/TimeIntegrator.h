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
#include "output/Output.h"
#include "output/sqlite/DatabaseConnection.h"

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
		bounds_(bounds), particles_(particles), conf_(conf), pusher_(pusher), tree_(tree)
	{}

	void addParticleOutput(output::Output<Vec>* o){
		per_particle_outputs_.push_back(o);
	}

	void addTimestepOutput(output::Output<Vec>* o){
		per_timestep_outputs_.push_back(o);
	}

	/**
	 * @brief Kick off the integration.
	 * @param precision		Precision to supply to the potential.
	 * @param output_every	Output every nth step.
	 * @param pos_file		File to output particle positions to (or NULL).
	 * @param vel_file		File to output particle velocities to (or NULL).
	 * @param energy_file	File to output energies to (or NULL).
	 */
	void start(potentials::Precision precision, unsigned int output_every,
			output::DatabaseConnection<Vec>& db){

		long int num_steps = conf_.getMaxTime() / conf_.getTimestep();
		for(long int i=0;i<num_steps;i++){
			std::pair<double, double> energies = pusher_.push_particles(particles_, tree_, bounds_, precision);
			bounds_.timestepOver();

			std::cout << "Timestep " << i << " of " << num_steps << " complete (" << ((float)i/(float)num_steps)*100 << "%)" << std::endl;
			if( (i%output_every) == 0){

				db.write_timestep_particles(i, particles_);
				db.writeEnergies(i, energies.first, energies.second);
			}
		}
	}

	/**
	 * @brief Destructor.
	 */
	~TimeIntegrator() {}

private:
	BoundaryConditions<Vec>& bounds_;
	std::vector<Particle<Vec>*>& particles_;
	const Configuration<Vec>& conf_;
	pusher::Pusher<Vec,Mat>& pusher_;
	Tree<Vec,Mat>& tree_;

	std::vector<output::Output<Vec>*> per_particle_outputs_;
	std::vector<output::Output<Vec>*> per_timestep_outputs_;
};

} /* namespace treecode */
#endif /* TIMEINTEGRATOR_H_ */
