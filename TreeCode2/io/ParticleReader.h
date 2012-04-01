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
	ParticleReader(const char *pos_file, const char *vel_file, double mass, int charge):
		pos_input_(NULL), vel_input_(NULL), mass_(mass), charge_(charge){

		if(pos_file == NULL)
			throw ReadError("Must specify a position file.");

		pos_input_ = new std::ifstream(pos_file);
		if(pos_input_->fail())
			throw ReadError("Failed to open position file.");

		if(vel_file != NULL){
			vel_input_ = new std::ifstream(vel_file);
			if(vel_input_->fail())
				throw ReadError("Failed to open velocity file.");
		}
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
	template<class OutputIterator>
	void readParticles(OutputIterator in, int timestep = 0){
		//Read at least timestep_ newlines. If we reach the end of the file,
		//throw an exception.
		int timesteps = 0;
		while(timesteps < timestep && !pos_input_->eof()){
			int c = pos_input_->get();
			if(c == '\n')
				timesteps ++;
		}
		if(pos_input_->eof())
			throw ReadError("Reached end of position file before reading the correct number of timesteps.");
		//Now do the same for the velocity file, if it exists
		if(vel_input_ != NULL){
			timesteps = 0;
			while(timesteps < timestep && !vel_input_->eof()){
				int c = vel_input_->get();
				if(c == '\n')
					timesteps ++;
			}
			if(pos_input_->eof())
				throw ReadError("Reached end of velocity file before reading the correct number of timesteps.");
		}

		Vec position = Vec::Zero();
		Vec velocity = Vec::Zero();
		int index = 0;	//When this reaches D, create a particle
		//Read D components, then create and add a particle
		while((*pos_input_) >> position[index]){

			//read in the velocity as well, if the file is specified
			if(vel_input_ != NULL){
				(*vel_input_) >> velocity[index];
				//If we have reached the end of the velocity file, throw
				//an exception, because it indicates there is the wrong
				//number of fields
				if(vel_input_->eof())
					throw ReadError("Reached end of file in velocities before positions.");
			}
			index = (index + 1) % D;
			if(index == 0){
				Particle<D> *p = new Particle<D>(charge_, mass_, position, velocity);
				//parts.push_back(p);
				in = p;
			}
			//If we have read the complete timestep, break.
			if(pos_input_->get() == '\n'){
				//We should also make sure that we are at the end of
				//the velocities file, to make sure they have
				//an equal number of fields.
				if(vel_input_ != NULL && vel_input_->get() != '\n')
					throw ReadError("Reached end of line in velocity file at different time to position file.");
				//We're done reading the line.
				break;
			}
		}
	}

private:
	std::ifstream *pos_input_, *vel_input_;
	double mass_;
	int charge_;
};

}
}

#endif /* PARTICLEREADER_H_ */
