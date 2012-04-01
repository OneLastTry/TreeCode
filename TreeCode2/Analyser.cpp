/*
 * Analyser.cpp
 *
 *  Created on: 1 Apr 2012
 *      Author: stefans
 */

#include <string>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <cmath>
#include <map>

#include <OptionParser.h>
#include <Particle.h>
#include <io/ParticleReader.h>

namespace po = boost::program_options;

typedef treecode::Particle<3> part_3d;
typedef std::vector<part_3d* >	part_list_3d;

int main(int argc, char **argv){
	using namespace treecode::io;

	OptionParser opts("Options");
	po::options_description actions("Possible actions");
	actions.add_options()
		("speeds", "Output speeds of all particles")
		("coll-time", "Produce a graph of average deviated angle against time")
		("radial-density", "Produce plot of radial density")
	;
	po::options_description compulsory("Compulsory arguments");
	compulsory.add_options()
		("pos-file,P", po::value<std::string>()->required(), "File containing position of particles")
		("vel-file,V", po::value<std::string>()->required(), "File containing velocities of particles")
		("dimensions,d", po::value<int>()->required(), "Dimensions in which we are working")
	;

	po::options_description optional("Optional (or conditional) arguments");
	optional.add_options()
		("charge,q", po::value<int>()->default_value(1), "Charge of particles")
		("mass,m", 	 po::value<double>()->default_value(1), "Mass of particles")
		("timestep,t", po::value<int>()->default_value(0), "Timestep to examine")
		("origin", po::value<Eigen::VectorXd>(), "Origin of system")
		("bin-width", po::value<double>(), "Bin width")
	;
	opts.add(actions).add(compulsory).add(optional);
	opts.parse(argc, argv);

	std::string pos_file = opts.get<std::string>("pos-file");
	std::string vel_file = opts.get<std::string>("vel-file");
	double mass = opts.get<double>("mass");
	int charge = opts.get<int>("charge");
	int timestep = opts.get<int>("timestep");

	if(opts.count("speeds")){
		ParticleReader<3> reader(pos_file.c_str(), vel_file.c_str(), mass, charge);
		part_list_3d parts = reader.readParticles(timestep);
		int size = parts.size();
		for(int i = 0; i < size; i++){
			std::cout << parts[i]->getVelocity().norm() << std::endl;
		}
	}else if(opts.count("coll-time")){
		ParticleReader<3> reader(pos_file.c_str(), vel_file.c_str(), mass, charge);
		//Read initial particles
		part_list_3d initial_parts = reader.readParticles(0);
		//Then just keep reading until we run out of particles
		part_list_3d parts;
		while(!reader.eof()){
			parts = reader.readParticles();
			//Average angle for this timestep
			double average = 0;
			int size = parts.size();
			for(int i = 0; i < size; i++){
				Eigen::Vector3d init_vel = initial_parts[i]->getVelocity().normalized();
				Eigen::Vector3d curr_vel = parts[i]->getVelocity().normalized();
				average += acos(init_vel.dot(curr_vel)) / size;
			}
			std::cout << average << std::endl;
		}
	}else if(opts.count("radial-density")){
		Eigen::VectorXd origin;
		double bin_width;
		try{
			 origin = opts.get<Eigen::VectorXd>("origin");
			 bin_width = opts.get<double>("bin-width");
		}catch(...){
			std::cerr << "You must specify the origin and bin width when plotting radial density" << std::endl;
			exit(1);
		}
		ParticleReader<3> reader(pos_file.c_str(), vel_file.c_str(), mass, charge);
		part_list_3d parts = reader.readParticles(timestep);
		//Bin index, density
		std::map<int, double> histogram;
		int size = parts.size();
		for(int i=0;i<size;i++){
			double r = (parts[i]->getPosition() - origin).norm();
			int bin = r / bin_width;
			//Volume of shell
			double vol = 4.0/3 * M_PI * ( pow(bin_width*(bin+1),3) - pow(bin_width*bin,3) );
			histogram[bin] += 1.0 / vol;
		}

		std::map<int,double>::iterator it;
		for ( it=histogram.begin() ; it != histogram.end(); it++ ){
		    std::cout << (it->first*bin_width)  << "\t" << it->second << std::endl;
		}
	}

	return 0;
}
