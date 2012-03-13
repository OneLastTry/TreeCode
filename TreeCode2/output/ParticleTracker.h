/*
 * ParticleTracker.h
 *
 *  Created on: 3 Mar 2012
 *      Author: stefans
 */

#ifndef PARTICLETRACKER_H_
#define PARTICLETRACKER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/foreach.hpp>
#include "../Particle.h"

namespace treecode{
namespace output{

template <class Vec>
class ParticleTracker{
public:
	enum token {
		POSITION,
		VELOCITY
	};

	ParticleTracker(std::string filename, const std::vector<Particle<Vec>*>& parts, token record):
		parts_(parts), record_(record){
		if(filename.compare("stdout") == 0)
			out_ = &std::cout;
		else if(filename.compare("stderr") == 0)
			out_ = &std::cerr;
		else
			out_ = new std::ofstream(filename.c_str());
	}

	virtual ~ParticleTracker(){
		if(out_->rdbuf() != std::cout.rdbuf() && out_->rdbuf() != std::cerr.rdbuf())
			delete out_;
	}


	virtual void output(){
		BOOST_FOREACH(Particle<Vec>* p, parts_){
			if(record_ == POSITION){
				for(unsigned int i = 0; i < p->getPosition().rows();i++)
					(*out_) << p->getPosition()[i] << "\t";
			}else if(record_ == VELOCITY){
				for(unsigned int i = 0; i < p->getVelocity().rows();i++)
					(*out_) << p->getVelocity()[i] << "\t";
			}
		}
		(*out_) << std::endl;
	}
protected:
	std::ostream* out_;
	const std::vector<Particle<Vec>* >& parts_;
	token record_;
};

}//output namespace
}//treecode namespace

#endif /* PARTICLETRACKER_H_ */
