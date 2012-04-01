/**
 * @file ParticleGenerator.cpp
 *
 * This source file is for a program to generate files containing
 * descriptions of particles. These files should then be read back
 * in and used as initial conditions.
 */

#include <boost/random.hpp>
#include <boost/program_options.hpp>

#include <distributions/dists.h>
#include <OptionParser.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

#include <distributions/Distribution.h>

namespace po = boost::program_options;


/**
 * Parse command line args, create distribution generators, and
 * output files
 */
int main(int argc, char **argv) {
	using namespace treecode::distribution;
	using boost::random::mt19937;
	using std::cout;

	//Variables to be set by program_options
	std::string pos_dist_str, vel_dist_str;
	std::string pos_out_str, vel_out_str;\
	Eigen::VectorXd origin_vector, const_vel_vector;
	double length;
	int num_particles;
	int random_seed;

	OptionParser opts("Options");
	opts.setUsesConfigFile(true);
	//Main variables. These don't depend on anything else, and are mandatory
	po::options_description main_opts("Main options");
	main_opts.add_options()
	    ("num-particles,N", po::value<int>(&num_particles)->required(), "Number of particles")
	    ("origin,O", 		po::value<Eigen::VectorXd>(&origin_vector)->required(), "Origin of system")
	    ("length,L", 		po::value<double>(&length)->required(), "Length of system")
	    ("pos-dist,p",		po::value<std::string>(&pos_dist_str)->required(), "Position distribution")
	    ("vel-dist,v",		po::value<std::string>(&vel_dist_str)->required(), "Velocity distribution")
	    ("pos-out,P",		po::value<std::string>(&pos_out_str)->required(), "Position output file")
	    ("vel-out,V",		po::value<std::string>(&vel_out_str)->required(), "Velocity output file")
	    ("seed,s",			po::value<int>(&random_seed)->required(), "Random seed")
	;

	//Options that only apply when using a maxwell distribution
	double temperature, mass;
	po::options_description maxwell_opts("Maxwell distribution options");
	maxwell_opts.add_options()
		("temperature,T", 	po::value<double>(&temperature), 	"Temperature")
		("mass,m",			po::value<double>(&mass), 	"Mass of particles")
	;

	//Options that only apply with a constant speed distribution
	po::options_description const_opts("Constant distribution options");
	const_opts.add_options()
		("velocity",	po::value<Eigen::VectorXd>(&const_vel_vector), "Constant velocity")
	;

	//Read options int vm
	opts.add(main_opts).add(maxwell_opts).add(const_opts);
	opts.parse(argc, argv);

	//Initialise output streams
	std::ofstream pos_out(pos_out_str.c_str());
	std::ofstream vel_out(vel_out_str.c_str());
	pos_out.precision(20);
	vel_out.precision(20);

	//Create distributions
	VectorDistribution* pos_dist;
	VectorDistribution* vel_dist;


	mt19937 rng(random_seed);
	//Parse the strings supplied, and generate the distributions
	if(pos_dist_str == "uniform")
		pos_dist = new UniformDistribution<mt19937>(rng, origin_vector, origin_vector.array() + length);
	else if(pos_dist_str == "spherical")
		pos_dist = new SphericalDistribution<mt19937>(rng, origin_vector.rows(), origin_vector, length);
	else{
		std::cerr << "Position distribution must be 'uniform' or 'spherical'." << std::endl;
		exit(1);
	}

	//Same, with velocity distributions.
	if(vel_dist_str == "maxwell"){
		if(!opts.count("mass") || !opts.count("temperature")){
			std::cerr << "Must specify temperature and mass for maxwell distribution" << std::endl;
			exit(1);
		}
		vel_dist = new MaxwellDistribution<mt19937>(rng, mass, temperature, origin_vector.rows());
	}
	else if(vel_dist_str == "constant"){
		if(!opts.count("velocity")){
			std::cerr << "Must specify a velocity with the constant velocity distribution" << std::endl;
			exit(1);
		}
		vel_dist = new ConstDistribution(const_vel_vector);
	}

	//Actually generate the particles
	for(int i=0;i<num_particles;i++){
		Eigen::VectorXd pos = pos_dist->getVector();
		Eigen::VectorXd vel = vel_dist->getVector();
		for(int j = 0;j<pos.rows();j++){
			pos_out << std::scientific << pos[j] << "\t";
			vel_out << std::scientific << vel[j] << "\t";
		}
	}
}
