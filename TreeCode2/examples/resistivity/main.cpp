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

#include <3d_typedefs.h>

#include <opt_parser/OptionParser.h>
#include <bounds/CullingBoundary.h>
#include <potentials/CoulombForceEField.h>
#include <macs/BarnesHutMAC.h>

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

void parse_cmd_line(int argc, char **argv,
		double &theta, double &fs, double &dt, double &max_time,
		double& length, unsigned int& num_parts, double& temperature){
	// Declare the supported options.
	OptionParser op(argc, argv);
	bool help;
	op <<
	    new BoolOption("--help" ,"-h", "Produce this message", help) <<
	    new ArgOption<double>("--length", 			"-l", 	"Length of one side of bounding box.", length) <<
	    new ArgOption<double>("--force-softening",	"-f",	"Force softening constant", fs) <<
	    new ArgOption<double>("--timestep", 		"-d",	"Timestep", dt) <<
	    new ArgOption<double>("--max-time", 		"-m", 	"Maximum time", max_time) <<
	    new ArgOption<double>("--theta", 			"-t", 	"Theta (MAC)", theta) <<
	    new ArgOption<unsigned int>("--number", 	"-n", 	"Number of each species", num_parts) <<
	    new ArgOption<double>("--temperature", 		"-T", 	"Temperature of plasma", temperature);

	unsigned int num_options = op.parse();

	if(help || num_options != op.size() - 1){
		op.display(std::cout);
		exit(0);
	}
}

int main(int argc, char **argv) {

	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;

	using boost::random::mt19937;
	mt19937 rng;

	double theta, force_softening, timestep, max_time, length, temperature;
	unsigned int num_particles;
	parse_cmd_line(argc, argv,
			theta, force_softening, timestep, max_time,
			length, num_particles, temperature);


	UniformDistribution3d 			position_dist(Vec(0,0,0), Vec(length, length, length));
	ConstDistribution3d				i_velocity_dist(Vec::Zero());
	MaxwellDistribution3d			e_velocity_dist(1,temperature);
	ConstantChargeDistribution3d	electron_charges(-1);
	ConstantChargeDistribution3d	ion_charges(1);

	int id = 0;
	vector<Particle3d*> parts;
	vector<Particle3d*> ions = Particle3d::generateParticles<mt19937>(num_particles, 1837, rng,
			position_dist, i_velocity_dist, ion_charges, id);
	vector<Particle3d*> electrons = Particle3d::generateParticles<mt19937>(num_particles, 1, rng,
				position_dist, e_velocity_dist, electron_charges, id);
	parts.insert(parts.end(), ions.begin(), ions.end());
	parts.insert(parts.end(), electrons.begin(), electrons.end());


	std::ofstream fout("count.csv");
	bool min_reset[] = {true, false, false};
	bool max_reset[] = {false, false, false};
	CullingBoundary<Vec,Mat,mt19937> bounds(Vec(0,0,0), length, e_velocity_dist, min_reset, max_reset, rng);

	CoulombForceEField<Vec, Mat> open_pot(force_softening, bounds, Vec(100, 0, 0));
	EwaldForce3d			periodic_pot(force_softening, bounds, 2.0 / length, 5, 5);
	InterpolatedEwaldSum3d	potential(force_softening, bounds, 20, periodic_pot, open_pot);
	potential.init();
	BarnesHutMAC<Vec,Mat>	mac(theta, bounds);
	LeapfrogPusher3d 		push(timestep, bounds, potential);
	Tree3d					tree(bounds, parts);
	TimeIntegrator3d		integrator(timestep, max_time, parts, tree, bounds, push, mac);

	integrator.setEnergyOutputFile("energies.csv");
	integrator.addParticleTracker(new ParticleTracker<Vec,Mat>("velocities.csv", electrons, ParticleTracker<Vec,Mat>::VELOCITY));

	push.init(parts, tree, quadrupole, mac);

	integrator.start(quadrupole, 1);

	return 0;
}
