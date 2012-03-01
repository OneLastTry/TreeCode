/*
 * main.cpp
 *
 *  Created on: 6 Dec 2011
 *      Author: stefans
 */

#include "Node.h"
#include "Particle.h"
#include "Configuration.h"
#include "bounds/OpenBoundary.h"
#include "bounds/PeriodicBoundary.h"
#include "Tree.h"
#include "distributions/UniformDistribution.h"
#include "distributions/MaxwellDistribution.h"
#include "distributions/ConstantChargeDistribution.h"
#include "distributions/ConstDistribution.h"
#include "distributions/SphericalDistribution.h"
#include "distributions/CylindricalDistribution.h"
#include "distributions/SinusoidalDistribution.h"
#include "pushers/LeapfrogPusher.h"
#include "potentials/CoulombForce.h"
#include "TimeIntegrator.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/chrono.hpp>
#include <Eigen/Dense>
#include <vector>

#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace treecode;
using namespace treecode::pusher;
using namespace treecode::potentials;
using namespace treecode::distribution;
using namespace std;

#ifdef __CDT_PARSER__
#define endl 0
#endif

//	high_resolution_clock::time_point start = high_resolution_clock::now();
//	milliseconds duration = boost::chrono::duration_cast<milliseconds>(high_resolution_clock::now() - start);
//	cout << "Tree build took " << duration << endl;

int main(void) {
	typedef Eigen::Vector3d ThreeVec;
	typedef Eigen::Matrix3d ThreeMat;

	using boost::chrono::time_point;
	using boost::chrono::high_resolution_clock;
	using boost::chrono::milliseconds;
	using boost::random::mt19937;

	using output::Output;


	int num_particles = 250;
	Configuration<ThreeVec> c(3, 0.3, 0.01, 100, 500, 0.1);

	mt19937 rng(2);
	ConstDistribution<mt19937, ThreeVec> electron_velocities(ThreeVec::Zero());
	ConstDistribution<mt19937, ThreeVec> proton_velocities(ThreeVec::Zero());

	UniformDistribution<mt19937,ThreeVec> particle_positions(ThreeVec::Zero(), ThreeVec(1,1,1));

	ConstantChargeDistribution<mt19937> proton_charges(+1);
	ConstantChargeDistribution<mt19937> electron_charges(-1);

	int id = 0;
	vector<Particle<ThreeVec>*> protons = Particle<ThreeVec>::generateParticles<mt19937>(num_particles, 1837, rng,
			particle_positions, proton_velocities, proton_charges, id);

	vector<Particle<ThreeVec>*> electrons = Particle<ThreeVec>::generateParticles<mt19937>(num_particles, 1, rng,
			particle_positions, electron_velocities, electron_charges, id);

	vector<Particle<ThreeVec>*> parts;
	parts.insert(parts.end(), protons.begin(), protons.end());
	parts.insert(parts.end(), electrons.begin(), electrons.end());

//	Output* energies = new Output("energies.csv");
//	energies->addRecord(Output::TIME);
//	energies->addRecord(Output::KINETIC_ENERGY);
//	energies->addRecord(Output::POTENTIAL_ENERGY);
//
//	Output* coordinates = new Output("coordinates.csv");
//	coordinates->addRecord(Output::COORDINATES);

	OpenBoundary<ThreeVec> ob(c);
	ob.init(parts);
	CoulombForceThreeD<ThreeVec, ThreeMat> coulomb_pot(c, ob);

//	PeriodicBoundary pb(c, Eigen::Vector3d(0,0,0), 1);
//	EwaldForce ewald_pot(c, pb, 2.0 / pb.getSize(), 5, 5);

	output::DatabaseConnection<ThreeVec> conn("test_db.db");
	conn.clear_database();
	conn.init_database(3);
	conn.write_sim_params(num_particles, c, ob);
	conn.write_init_particles(parts);

	cout << "Generating Ewald grid" << endl;
//	InterpolatedEwaldSum interpolated(c, pb, 20, ewald_pot, coulomb_pot);
//	interpolated.init();

	cout << "Performing initial tree build" << endl;
	Tree<ThreeVec, ThreeMat> tr(c, ob, parts);
	tr.rebuild();

	cout << "Initialising leapfrog pusher" << endl;
	LeapfrogPusher<ThreeVec, ThreeMat> leapfrog(c, ob, coulomb_pot);
	leapfrog.init(parts, tr, Precision::quadrupole);

	cout << "Beginning simulation" << endl;
	TimeIntegrator<ThreeVec, ThreeMat> time_integrator(c, parts, tr, ob, leapfrog);
	time_integrator.start(Precision::quadrupole, 10, conn);

	Particle<ThreeVec>::deleteParticles(parts);

//	delete energies;
//	delete coordinates;

	return 0;
}
