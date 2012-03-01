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

#include <bounds/CountingPeriodicBounds.h>
#include <potentials/CoulombForceEField.h>

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

#define E0 8.85E-12
#define E_CHARGE 1.62E-19
#define M_e 9.11E-31

Configuration3d parse_cmd_line(int argc, char **argv, double& length, unsigned int& num_parts, std::string& dbname){
	// Declare the supported options.
	po::options_description desc("Allowed options");
	double fs, dt, max_time, theta, param, temperature;
	desc.add_options()
	    ("help,h", 				"Produce this message")
	    ("length,l", 		 	po::value<double>(&length), "Length of one side of bounding box.")
	    ("force-softening,f",	po::value<double>(&fs), "Force softening constant")
	    ("timestep,d", 		po::value<double>(&dt), "Timestep")
	    ("max-time,m", 			po::value<double>(&max_time), "Maximum time")
	    ("theta,t", 			po::value<double>(&theta), "Theta (MAC)")
	    ("param,p", 			po::value<double>(&param), "Plasma parameter (bigger implies more ideal)")
	    ("number,n", 			po::value<unsigned int>(&num_parts), "Number of each species")
	    ("database,b", 			po::value<std::string>(&dbname), "Database file name")
	    ("temperature,T", 		po::value<double>(&temperature), "Temperature of plasma")
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

	double length;
	unsigned int num_particles;
	std::string dbname;
	Configuration3d c = parse_cmd_line(argc, argv, length, num_particles, dbname);

	UniformDistribution3d 			position_dist(Vec(0,0,0), Vec(length, length, length));
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


	double density = 1E14; //cm3
	double T = 100;	//eV
	double Nd = c.getDensity(); //Plasma param
	double E_field = 1E2; // V/m

	double debye_length = sqrt(E0 * T / (density * 1E6 * E_CHARGE));
	double thermal_velocity = sqrt(T * E_CHARGE / M_e);
	double plasma_freq = sqrt(density * 1E6 * E_CHARGE * E_CHARGE / (M_e*E0));
	double E_field_prime = E_field * E_CHARGE  / (M_e * debye_length * plasma_freq * plasma_freq);

	cout << "Reference plasma has: n_e = " << density << " / cm3, T = " <<
			T << "eV, Nd = " << Nd << endl;
	cout << "Therefore, the Debye length is: " << debye_length << "m" << endl;
	cout << "The thermal velocity is: " << thermal_velocity << "m/s" << endl;
	cout << "The plasma frequency is: " << plasma_freq << "s-1" << endl;
	cout << "An electric field of " << E_field << "V/m is therefore " <<
			E_field_prime << " in the dimensionless units" << endl;


	std::ofstream fout("test");
	CountingPeriodicBounds<Vec> bounds(c, Vec(0,0,0), length, fout);

	CoulombForceEField<Vec, Mat> open_pot(c, bounds, Vec(E_field_prime, 0, 0));
	EwaldForce3d			periodic_pot(c, bounds, 2.0 / length, 5, 5);
	InterpolatedEwaldSum3d	potential(c, bounds, 20, periodic_pot, open_pot);
	LeapfrogPusher3d 		push(c, bounds, potential);
	Tree3d					tree(c, bounds, parts);
	DatabaseConnection3d	db(dbname);
	TimeIntegrator3d		integrator(c, parts, tree, bounds, push);

	db.clear_database();
	db.init_database(3);
	db.write_sim_params(num_particles, c, bounds);
	db.write_init_particles(parts);

	push.init(parts, tree, quadrupole);

	integrator.start(quadrupole, 10000, db);

	return 0;
}
