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



Configuration3d parse_cmd_line(int argc, char **argv, double& radius, unsigned int& num_parts, std::string& dbname){
	// Declare the supported options.
	po::options_description desc("Allowed options");
	double fs, dt, max_time, theta, param;
	desc.add_options()
	    ("help,h", 				"Produce this message")
	    ("radius,r", 		 	po::value<double>(&radius), "Radius of starting plasma")
	    ("force-softening,f",	po::value<double>(&fs), "Force softening constant")
	    ("timestep,d", 		po::value<double>(&dt), "Timestep")
	    ("max-time,m", 			po::value<double>(&max_time), "Maximum time")
	    ("theta,t", 			po::value<double>(&theta), "Theta (MAC)")
	    ("param,p", 			po::value<double>(&param), "Plasma parameter (bigger implies more ideal)")
	    ("number,n", 			po::value<unsigned int>(&num_parts), "Number of each species")
	    ("database,b", 			po::value<std::string>(&dbname), "Database file name")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help") || vm.size() != 8) {
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

	double radius;
	unsigned int num_particles;
	std::string dbname;
	Configuration3d c = parse_cmd_line(argc, argv, radius, num_particles, dbname);

	SpherialDistribution3d			position_dist(3, Vec::Zero(), radius);
	ConstDistribution3d				i_velocity_dist(Vec::Zero());
	MaxwellDistribution3d			e_velocity_dist(1,1);
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

	OpenBoundary3d 			bounds(c);
	CoulombForce3d 			potential(c, bounds);
	LeapfrogPusher3d 		push(c, bounds, potential);
	Tree3d					tree(c, bounds, parts);
	DatabaseConnection3d	db(dbname);
	TimeIntegrator3d		integrator(c, parts, tree, bounds, push);

	db.clear_database();
	db.init_database(3);
	db.write_sim_params(num_particles, c, bounds);
	db.write_init_particles(parts);

	bounds.init(parts);
	push.init(parts, tree, quadrupole);

	integrator.start(quadrupole, 10, db);


	cout << radius << endl;

	return 0;
}
