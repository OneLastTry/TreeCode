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

#include <opt_parser/OptionParser.h>
#include <Configuration.h>
#include <Particle.h>
#include <Node.h>
#include <Tree.h>
#include <TimeIntegrator.h>
#include <macs/BarnesHutMAC.h>
#include <output/ParticleTracker.h>

#include <Eigen/Dense>
#include <vector>
#include <string>

using namespace std;
using namespace treecode;

namespace po = boost::program_options;



Configuration3d parse_cmd_line(int argc, char **argv, double& length, unsigned int& num_parts, double& proportion, double& wavelengths, double& temperature){

	OptionParser op(argc, argv);
	bool help;
	double timestep, theta, max_time, param, fs;

	op << new ArgOption<double>("--length", "-l", "Length of system (in debye lengths)", length) <<
			new ArgOption<double>("--timestep", "-dt", "Individual timestep.", timestep) <<
			new ArgOption<double>("--softening", "-fs", "Force softening parameter.", fs) <<
			new ArgOption<double>("--max-time", "-mt", "Maximum time to run to.", max_time) <<
			new ArgOption<double>("--theta", "-t", "Critical opening angle.", theta) <<
			new ArgOption<double>("--param", "-p", "Plasma parameter (larger implies more ideal)", param) <<
			new ArgOption<unsigned int>("--number", "-n", "Number of each species.", num_parts) <<
			new ArgOption<double>("--proportion", "-prop", "Proportion of particles in sinusoidal perturbation (0.0-1.0).", proportion) <<
			new ArgOption<double>("--wavelengths", "-w", "Number of wavelengths in perturbation.", wavelengths) <<
			new ArgOption<double>("--temperature", "-temp", "Temperature of plasma.", temperature) <<
			new BoolOption("--help", "-h", "This help text.", help);
	unsigned int num_options = op.parse();

	if(help || num_options != op.size() - 1){
		op.display();
		exit(0);
	}

	return Configuration3d(3, theta, timestep, max_time, param, fs);

}

int main(int argc, char **argv) {
	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;

	using boost::mt19937;
	mt19937 rng;


	double length, wavelengths, proportion, temperature;
	unsigned int num_particles;
	std::string dbname;
	Configuration3d c = parse_cmd_line(argc, argv, length, num_particles, proportion, wavelengths, temperature);

	Vec origin(-length/2,-length/2,-length/2);
	Vec max(length/2, length/2, length/2);
	UniformDistribution3d 			position_dist(origin, max);
	SinusoidalDistribution3d		perturbed_dist(0, origin, max, wavelengths, M_PI_2);
	ConstDistribution3d				i_velocity_dist(Vec::Zero());
	MaxwellDistribution3d			e_velocity_dist(1, temperature);
	ConstantChargeDistribution3d	electron_charges(-1);
	ConstantChargeDistribution3d	ion_charges(1);

	int id = 0;
	vector<Particle3d*> parts;
	vector<Particle3d*> ions = Particle3d::generateParticles<mt19937>(num_particles, 100000, rng,
			position_dist, i_velocity_dist, ion_charges, id);
	vector<Particle3d*> electrons = Particle3d::generateParticles<mt19937>(num_particles * (1 - proportion), 1, rng,
				position_dist, e_velocity_dist, electron_charges, id);
	vector<Particle3d*> perturbed_electrons = Particle3d::generateParticles<mt19937>(num_particles*proportion, 1, rng,
					perturbed_dist, e_velocity_dist, electron_charges, id);

	parts.insert(parts.end(), ions.begin(), ions.end());
	parts.insert(parts.end(), electrons.begin(), electrons.end());
	parts.insert(parts.end(), perturbed_electrons.begin(), perturbed_electrons.end());

	PeriodicBoundary3d		bounds(c, origin, length);
	BarnesHutMAC<Vec,Mat>		mac(0.7, bounds);
	CoulombForce3d 			open_pot(c, bounds);
	EwaldForce3d			periodic_pot(c, bounds, 2.0 / length, 10, 10);
	InterpolatedEwaldSum3d	potential(c, bounds, 25, periodic_pot, open_pot);
	potential.init();
	LeapfrogPusher3d 		push(c, bounds, potential);
	Tree3d					tree(c, bounds, parts);
	TimeIntegrator3d		integrator(c, parts, tree, bounds, push, mac);
	integrator.setEnergyOutputFile("energies.csv");
	integrator.addParticleTracker(new ParticleTracker<Vec,Mat>("positions.csv", parts, ParticleTracker<Vec,Mat>::POSITION));
	integrator.addParticleTracker(new ParticleTracker<Vec,Mat>("velocities.csv", parts, ParticleTracker<Vec,Mat>::VELOCITY));

	cout << "Initialising pusher" << endl;

	push.init(parts, tree, quadrupole, mac);

	integrator.start(quadrupole, 1);

	return 0;
}
