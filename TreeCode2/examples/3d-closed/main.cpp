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



Configuration3d parse_cmd_line(int argc, char **argv, double& length, unsigned int& num_parts, std::string& dbname){
	// Declare the supported options.
	po::options_description desc("Allowed options");
	double fs, dt, max_time, theta, param;
	desc.add_options()
	    ("help,h", 				"Produce this message")
	    ("length,l", 		 	po::value<double>(&length), "Length of one side of bounding box.")
	    ("force-softening,f",	po::value<double>(&fs), "Force softening constant")
	    ("timestep,d", 		po::value<double>(&dt), "Timestep")
	    ("max-time,m", 			po::value<double>(&max_time), "Maximum time")
	    ("theta,t", 			po::value<double>(&theta), "Theta (MAC)")
	    ("param,p", 			po::value<double>(&param), "Plasma parameter (bigger implies more ideal)")
	    ("number,n", 			po::value<unsigned int>(&num_parts), "Number of each species")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help") || vm.size() != desc.options().size() - 1) {
	    cerr << desc << endl;
	    cerr << "All options must be specified!" << endl;
	    exit(1);
	}

	return Configuration3d(3, theta, dt, max_time, param, fs);
}

int main(int argc, char **argv) {
	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;

	using boost::random::mt19937;
	mt19937 rng;

	double length;
	unsigned int num_particles;
	std::string dbname;
	Configuration3d c = parse_cmd_line(argc, argv, length, num_particles, dbname);

	UniformDistribution3d 			position_dist(Vec(0,0,0), Vec(length, length, length));
	ConstDistribution3d				i_velocity_dist(Vec::Zero());
	ConstDistribution3d				e_velocity_dist(Vec::Zero());
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

	PeriodicBoundary3d		bounds(c, Vec(0,0,0), length);
	CoulombForce3d 			open_pot(c, bounds);
	EwaldForce3d			periodic_pot(c, bounds, 2.0 / length, 5, 5);
	InterpolatedEwaldSum3d	potential(c, bounds, 50, periodic_pot, open_pot);
	LeapfrogPusher3d 		push(c, bounds, potential);
	Tree3d					tree(c, bounds, parts);
	TimeIntegrator3d		integrator(c, parts, tree, bounds, push);

	push.init(parts, tree, quadrupole);

	integrator.start(quadrupole, 5);

	return 0;
}
