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
#include <Configuration.h>
#include <Particle.h>
#include <Node.h>
#include <Tree.h>
#include <TimeIntegrator.h>

#include <distributions/ParticleAvoidingDistribution.h>
#include <distributions/UniformDistribution.h>
#include <distributions/MaxwellDistribution.h>
#include <distributions/ConstantChargeDistribution.h>

#include <output/ParticleTracker.h>

#include <bounds/PeriodicBoundary.h>
#include <bounds/CountingPeriodicBounds.h>

#include <drude/BigParticle.h>
#include <drude/DrudeMAC.h>
#include <drude/CullingDrudePusher.h>

#include <Eigen/Dense>
#include <vector>
#include <string>

using namespace std;
using namespace treecode;

namespace po = boost::program_options;

typedef Eigen::Vector2d Vec;
typedef Eigen::Matrix2d Mat;

Configuration<Vec> parse_cmd_line(int argc, char **argv, unsigned int& num_parts){
//	OptionParser op(argc, argv);
//	double dt, max_time, theta;
//	bool display_help;
//			op << new ArgOption<double>("--timestep", "-dt", "Individual timestep.", dt) <<
//			new ArgOption<double>("--max-time", "-mt", "Maximum time to run to.", max_time) <<
//			new ArgOption<double>("--theta", "-t", "Critical opening angle.", theta) <<
//			new ArgOption<unsigned int>("--number", "-n", "Number of each species.", num_parts) <<
//			new BoolOption("--help", "-h", "This help text.", display_help);
//	unsigned int num_opts = op.parseOptions();
//
//	//If the user specified -h or not all options are set, display usage.
//	if(display_help || num_opts != op.numOptions() - 1){
//		op.display(std::cout);
//		exit(0);
//	}
	num_parts = 10;
	return Configuration<Vec>(2,0.0, 0.01, 100, 0, 0);
}

int main(int argc, char **argv) {
	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;

	typedef Particle<Vec> Particle_t;
	typedef BigParticle<Vec> BigParticle_t;


	using boost::mt19937;
	mt19937 rng;


	double temperature = 1;
	unsigned int num_particles;
	std::string dbname;
	Configuration<Vec> c = parse_cmd_line(argc, argv, num_particles);


	double length = 100;

	Vec origin(-length/2,-length/2);
	Vec max(length/2, length/2);
	MaxwellDistribution<mt19937,Vec>		e_velocity_dist(1, temperature);
	MaxwellDistribution<mt19937,Vec>		i_velocity_dist(1837000, temperature);
	ConstantChargeDistribution<mt19937>		electron_charges(-1);
	ConstantChargeDistribution<mt19937>		ion_charges(1);

	int id = 0;
	vector<Particle_t*> parts;
	vector<BigParticle_t*> big_ions;
	vector<Particle_t*> electrons;

	for(double x = -45; x < 50; x+= 15){
		for(double y = -45; y < 50; y+= 15){
			big_ions.push_back(new BigParticle<Vec>(+1, 100000, Vec(x,y), Vec(0,0), id, 5));
		}
	}

	ParticleAvoidingDistribution<mt19937, Vec, BigParticle_t> pos_vec(origin, max, big_ions);
	electrons = Particle<Vec>::generateParticles(num_particles, 1, rng, pos_vec, e_velocity_dist, electron_charges, id);


	parts.insert(parts.end(), big_ions.begin(), big_ions.end());
	parts.insert(parts.end(), electrons.begin(), electrons.end());

	CountingPeriodicBounds<Vec>	bounds(c, origin, length, cerr);
	CullingDrudePusher<Vec, Mat, mt19937> 	push(c, bounds, e_velocity_dist, rng);
	DrudeMAC<Vec, Mat>		mac(c.getTimestep());
	Tree<Vec,Mat>			tree(c, bounds, parts);
	TimeIntegrator<Vec,Mat>	integrator(c, electrons, tree, bounds, push, mac);
	integrator.setEnergyOutputFile("energies.csv");


	integrator.addParticleTracker(new ParticleTracker<Vec>("positions.csv", parts, ParticleTracker<Vec>::POSITION));
	integrator.addParticleTracker(new ParticleTracker<Vec>("velocities.csv", parts, ParticleTracker<Vec>::VELOCITY));


	integrator.start(quadrupole, 1);

	return 0;
}
