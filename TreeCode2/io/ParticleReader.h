/*
 * ParticleReader.h
 *
 *  Created on: 28 Mar 2012
 *      Author: stefans
 */

#ifndef PARTICLEREADER_H_
#define PARTICLEREADER_H_

#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <Eigen/Dense>
#include "../Particle.h"

namespace treecode{
namespace io{

class ReadError : public std::runtime_error{
public:
	ReadError(const char *msg):std::runtime_error(msg){}
};

template <int D>
class ParticleReader{
	typedef Eigen::Matrix<double, D, 1> Vec;

public:
	/**
	 * @brief Instantiate a new ParticleReader object, for particle input.
	 *
	 * If vel_file is not supplied, the velocities are initialised to
	 * a zero vector.
	 *
	 * @throw ReadError Thrown when both @p pos_file is NULL.
	 * @param pos_file	Location of position input file.
	 * @param vel_file	Location of velocity input file.
	 * @param mass		Mass of all particles to be instantiated.
	 * @param charge	Charge of all particles to be instantiated.
	 * @param timestep  Timestep to read.
	 */
	ParticleReader(std::string pos_file, std::string vel_file, double mass, int charge):
		pos_input_(NULL), vel_input_(NULL), mass_(mass), charge_(charge){

		pos_input_ = new std::ifstream(pos_file.c_str());
		vel_input_ = new std::ifstream(vel_file.c_str());

		if(pos_input_->fail())
			throw ReadError("Failed to open position file.");
		if(vel_input_->fail())
			throw ReadError("Failed to open velocity file.");
	}

	~ParticleReader(){
		delete pos_input_;
		delete vel_input_;
	}

	/**
	 * @brief Read particles from position and velocity files.
	 *
	 * The files themselves simply have the components of the
	 * position or velocity vector, separated by tabs, with
	 * each timestep on a different line.
	 *
	 * If the files don't have enough timesteps, a ReadError is
	 * thrown. If the files have differing number of fields, then
	 * a ReadError is thrown.
	 *
	 * @param timestep Timestep to examine
	 * @return Vector of Particle%s.
	 */
	std::vector<Particle<D>* > readParticles(int timestep_offset = 0){
		std::vector<Particle<D>* > parts;

		//Skip timestep_offset lines
		while(timestep_offset-- > 0 && !pos_input_->eof() && !vel_input_->eof()){
			pos_input_->ignore(10000000, '\n');
			pos_input_->unget();
			vel_input_->ignore(10000000, '\n');
			vel_input_->unget();
			if(pos_input_->get() == '\n' &&  vel_input_->get() == '\n')
				timestep_offset--;
		}

		//Read entire line into string
		std::string pos_line, vel_line;
		std::getline(*pos_input_, pos_line);
		std::getline(*vel_input_, vel_line);
		//Create stringstreams
		std::stringstream pos_ss(pos_line);
		std::stringstream vel_ss(vel_line);

		//Read into vectors, create new particle
		Vec pos, vel;
		int index = 0;
		while(pos_ss >> pos[index] && vel_ss >> vel[index]){
			index = (index + 1) % D;
			if(index == 0)
				parts.push_back(new treecode::Particle<D>(charge_, mass_, pos, vel));
		}
		return parts;
	}

	bool eof(){
		return pos_input_->eof() || vel_input_->eof();
	}

private:
	std::ifstream *pos_input_, *vel_input_;
	double mass_;
	int charge_;
};

}
}

#endif /* PARTICLEREADER_H_ */
