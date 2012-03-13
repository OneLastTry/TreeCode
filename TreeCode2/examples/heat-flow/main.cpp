/**
 * @file
 * @brief A sample program using open boundaries in 3D.
 *
 *
 */

#include <cstdlib>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/foreach.hpp>

#include <3d_typedefs.h>

#include <potentials/CoulombForceEField.h>
#include <bounds/CountingPeriodicBounds.h>
#include <opt_parser/OptionParser.h>
#include <Configuration.h>
#include <Particle.h>
#include <Node.h>
#include <Tree.h>
#include <TimeIntegrator.h>

#include <Eigen/Dense>
#include <vector>
#include <string>

using namespace std;
using namespace treecode;

namespace po = boost::program_options;



Configuration3d parse_cmd_line(int argc, char **argv, double& length, unsigned int& num_parts, double& temperature){
	OptionParser op(argc, argv);
	double dt, max_time, theta, plasma_param, force_softening;
	bool display_help;
	op << new ArgOption<double>("--length", "-l", "Length of system (in debye lengths)", length) <<
			new ArgOption<double>("--timestep", "-dt", "Individual timestep.", dt) <<
			new ArgOption<double>("--softening", "-fs", "Force softening parameter.", force_softening) <<
			new ArgOption<double>("--max-time", "-mt", "Maximum time to run to.", max_time) <<
			new ArgOption<double>("--theta", "-t", "Critical opening angle.", theta) <<
			new ArgOption<double>("--param", "-p", "Plasma parameter (larger implies more ideal)", plasma_param) <<
			new ArgOption<unsigned int>("--number", "-n", "Number of each species.", num_parts) <<
			new ArgOption<double>("--temperature", "-temp", "Temperature of plasma.", temperature) <<
			new BoolOption("--help", "-h", "This help text.", display_help);
	unsigned int num_opts = op.parseOptions();

	//If the user specified -h or not all options are set, display usage.
	if(display_help || num_opts != op.numOptions() - 1){
		op.display(std::cout);
		exit(0);
	}

	return Configuration3d(3,theta, dt, max_time, plasma_param, force_softening);
}

int main(int argc, char **argv) {
	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;

	using boost::mt19937;
	mt19937 rng;


	double length, temperature;
	unsigned int num_particles;
	Configuration3d c = parse_cmd_line(argc, argv, length, num_particles, temperature);


	Vec origin(-length/2,-length/2,-length/2);
	Vec max(length/2, length/2, length/2);
	UniformDistribution3d 			position_dist(origin, max);
	MaxwellDistribution3d			e_velocity_dist(1, temperature);
	MaxwellDistribution3d			i_velocity_dist(1837, temperature);
	ConstantChargeDistribution3d	electron_charges(-1);
	ConstantChargeDistribution3d	ion_charges(1);

	int id = 0;
	vector<Particle3d*> parts;
	vector<Particle3d*> ions = Particle3d::generateParticles<mt19937>(
			num_particles,
			1837,
			rng,
			position_dist,
			i_velocity_dist,
			ion_charges,
			id);
	vector<Particle3d*> electrons = Particle3d::generateParticles<mt19937>(
			num_particles,
			1,
			rng,
			position_dist,
			e_velocity_dist,
			electron_charges,
			id);

	parts.insert(parts.end(), ions.begin(), ions.end());
	parts.insert(parts.end(), electrons.begin(), electrons.end());

	MaxwellDistribution3d	perturbed_velocity_dist(1, 1);
	BOOST_FOREACH(Particle3d* p, electrons){
		if(p->getPosition().norm() < length/10)
			p->setVelocity(perturbed_velocity_dist.getVector(rng));
	}

	PeriodicBoundary3d		bounds(c, origin, length);
	CoulombForce3d			open_pot(c, bounds);
	EwaldForce3d			periodic_pot(c, bounds, 2.0 / length, 10, 10);
	InterpolatedEwaldSum3d	potential(c, bounds, 1, periodic_pot, open_pot);
	potential.init();
	LeapfrogPusher3d 		push(c, bounds, potential);
	Tree3d					tree(c, bounds, parts);
	TimeIntegrator3d		integrator(c, parts, tree, bounds, push);
	integrator.setEnergyOutputFile("energies.csv");

	integrator.addParticleTracker(new ParticleTracker<Vec>("positions.csv", parts, ParticleTracker<Vec>::POSITION));
	integrator.addParticleTracker(new ParticleTracker<Vec>("velocities.csv", parts, ParticleTracker<Vec>::VELOCITY));

	cout << "Initialising pusher" << endl;

	push.init(parts, tree, quadrupole);

	integrator.start(quadrupole, 1);

	return 0;
}
