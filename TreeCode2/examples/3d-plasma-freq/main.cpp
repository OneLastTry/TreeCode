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

#include <opt_parser/OptionParser.h>
#include <Particle.h>
#include <Node.h>
#include <Tree.h>
#include <TimeIntegrator.h>
#include <macs/BarnesHutMAC.h>
#include <output/CoordTracker.h>

#include <bounds/PeriodicBoundary.h>
#include <bounds/OpenBoundary.h>

#include <potentials/CoulombForce.h>
#include <potentials/EwaldForce.h>
#include <potentials/InterpolatedEwaldSum.h>
#include <pushers/LeapfrogPusher.h>

#include <distributions/UniformDistribution.h>
#include <distributions/MaxwellDistribution.h>
#include <distributions/SinusoidalDistribution.h>
#include <distributions/ConstDistribution.h>
#include <distributions/ConstantChargeDistribution.h>
#include <distributions/SphericalDistribution.h>

#include <Eigen/Dense>
#include <vector>
#include <string>

using namespace std;
using namespace treecode;

namespace po = boost::program_options;

typedef Eigen::Matrix<double, 3, 1> Vec;

void parse_cmd_line(int argc, char **argv,
		double& timestep, double &theta, double &max_time, double &fs,
		double& length, unsigned int& num_parts, double& proportion, double& wavelengths, double& temperature, bool& is_periodic){

	OptionParser op(argc, argv);
	bool help = false, periodic = false, open = false;

	op << new ArgOption<double>("--length", "-l", "Length of system (in debye lengths)", length) <<
			new ArgOption<double>("--timestep", "-dt", "Individual timestep.", timestep) <<
			new ArgOption<double>("--softening", "-fs", "Force softening parameter.", fs) <<
			new ArgOption<double>("--max-time", "-mt", "Maximum time to run to.", max_time) <<
			new ArgOption<double>("--theta", "-t", "Critical opening angle.", theta) <<
			new ArgOption<unsigned int>("--number", "-n", "Number of each species.", num_parts) <<
			new ArgOption<double>("--proportion", "-prop", "Proportion of particles in sinusoidal perturbation (0.0-1.0).", proportion) <<
			new ArgOption<double>("--wavelengths", "-w", "Number of wavelengths in perturbation.", wavelengths) <<
			new ArgOption<double>("--temperature", "-temp", "Temperature of plasma.", temperature) <<
			new BoolOption("--open", "-o", "Use open boundary conditions", open) <<
			new BoolOption("--periodic", "-p", "Use periodic boundary conditions", periodic) <<
			new BoolOption("--help", "-h", "This help text.", help);
	unsigned int num_options = op.parse();

	if(!periodic && !open){
		op.display(std::cerr);
		exit(1);
	}else if(periodic && open){
		std::cerr << "Must specify only one of --periodic or --open." << std::endl;
		exit(1);
	}
	is_periodic = periodic;

	if(help || num_options != op.size() - 2){
		op.display();
		exit(0);
	}
}

void periodic_start(double length, double wavelengths, double proportion,
		double temperature, double timestep, double theta, double max_time,
		double force_softening, unsigned int num_particles){
	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;

	using boost::mt19937;
	mt19937 rng;


	Vec origin(-length/2,-length/2,-length/2);
	Vec max(length/2, length/2, length/2);
	UniformDistribution<mt19937, 3> 			position_dist(origin, max);
	SinusoidalDistribution<mt19937, 3>		perturbed_dist(0, origin, max, wavelengths, M_PI_2);
	ConstDistribution<mt19937, 3>				i_velocity_dist(Vec::Zero());
	MaxwellDistribution<mt19937, 3>			e_velocity_dist(1, temperature);
	ConstantChargeDistribution<mt19937>	electron_charges(-1);
	ConstantChargeDistribution<mt19937>	ion_charges(1);

	int id = 0;
	vector<Particle<3>*> parts;
	vector<Particle<3>*> ions = Particle<3>::generateParticles<mt19937>(num_particles, 100000, rng,
			position_dist, i_velocity_dist, ion_charges, id);
	vector<Particle<3>*> electrons = Particle<3>::generateParticles<mt19937>(num_particles * (1 - proportion), 1, rng,
				position_dist, e_velocity_dist, electron_charges, id);
	vector<Particle<3>*> perturbed_electrons = Particle<3>::generateParticles<mt19937>(num_particles*proportion, 1, rng,
					perturbed_dist, e_velocity_dist, electron_charges, id);

	parts.insert(parts.end(), ions.begin(), ions.end());
	parts.insert(parts.end(), electrons.begin(), electrons.end());
	parts.insert(parts.end(), perturbed_electrons.begin(), perturbed_electrons.end());

	PeriodicBoundary<3>		bounds(origin, length);
	BarnesHutMAC<3>		mac(theta, bounds);
	CoulombForceThreeD<3> 			open_pot(force_softening, bounds);
	EwaldForce<3>			periodic_pot(force_softening, bounds, 2.0 / length, 5, 5);
	InterpolatedEwaldSum<3>	potential(force_softening, bounds, 20, periodic_pot, open_pot);
	potential.init();
	LeapfrogPusher<3> 		push(timestep, bounds, potential);
	Tree<3>					tree(bounds, parts);
	TimeIntegrator<3>		integrator(timestep, max_time, parts, tree, bounds, push, mac);
	integrator.setEnergyOutputFile("energies.csv");
	integrator.addParticleTracker(new CoordTracker<3>("positions.csv", parts, CoordTracker<3>::POSITION));
	integrator.addParticleTracker(new CoordTracker<3>("velocities.csv", parts, CoordTracker<3>::VELOCITY));

	cout << "Initialising pusher" << endl;

	push.init(parts, tree, quadrupole, mac);

	integrator.start(quadrupole, 1);
}

void open_start(double length, double wavelengths, double proportion,
		double temperature, double timestep, double theta, double max_time,
		double force_softening, unsigned int num_particles){
	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;

	using boost::mt19937;
	mt19937 rng;


	Vec origin = Vec::Zero();
	SphericalDistribution<mt19937, 3>			position_dist(3, origin, length);
	SphericalDistribution<mt19937, 3>			perturbed_dist(3, origin, length/20);
	ConstDistribution<mt19937, 3>				i_velocity_dist(Vec::Zero());
	MaxwellDistribution<mt19937, 3>			e_velocity_dist(1, temperature);
	ConstantChargeDistribution<mt19937>	electron_charges(-1);
	ConstantChargeDistribution<mt19937>	ion_charges(1);

	int id = 0;
	vector<Particle<3>*> parts;
	vector<Particle<3>*> ions = Particle<3>::generateParticles<mt19937>(num_particles, 100000, rng,
			position_dist, i_velocity_dist, ion_charges, id);
	vector<Particle<3>*> electrons = Particle<3>::generateParticles<mt19937>(num_particles * (1 - proportion), 1, rng,
				position_dist, e_velocity_dist, electron_charges, id);
	vector<Particle<3>*> perturbed_electrons = Particle<3>::generateParticles<mt19937>(num_particles*proportion, 1, rng,
					perturbed_dist, e_velocity_dist, electron_charges, id);

	parts.insert(parts.end(), ions.begin(), ions.end());
	parts.insert(parts.end(), electrons.begin(), electrons.end());
	parts.insert(parts.end(), perturbed_electrons.begin(), perturbed_electrons.end());

	OpenBoundary<3>			bounds;
	bounds.init(parts);
	BarnesHutMAC<3>	mac(theta, bounds);
	CoulombForceThreeD<3> 			potential(force_softening, bounds);
	LeapfrogPusher<3> 		push(timestep, bounds, potential);
	Tree<3>					tree(bounds, parts);
	TimeIntegrator<3>		integrator(timestep, max_time, parts, tree, bounds, push, mac);
	integrator.setEnergyOutputFile("energies.csv");
	integrator.addParticleTracker(new CoordTracker<3>("positions.csv", parts, CoordTracker<3>::POSITION));
	integrator.addParticleTracker(new CoordTracker<3>("velocities.csv", parts, CoordTracker<3>::VELOCITY));

	cout << "Initialising pusher" << endl;

	push.init(parts, tree, quadrupole, mac);

	integrator.start(quadrupole, 1);
}

int main(int argc, char **argv) {

	double length, wavelengths, proportion, temperature, timestep, theta, max_time, force_softening;
	unsigned int num_particles;
	bool periodic;
	parse_cmd_line(argc, argv, timestep, theta, max_time, force_softening,
			length, num_particles, proportion, wavelengths, temperature, periodic);

	if(periodic)
		periodic_start(length, wavelengths, proportion, temperature, timestep, theta, max_time, force_softening, num_particles);
	else
		open_start(length, wavelengths, proportion, temperature, timestep, theta, max_time, force_softening, num_particles);

	return 0;
}
