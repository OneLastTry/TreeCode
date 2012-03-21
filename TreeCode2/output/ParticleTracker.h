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

template <class Vec, class Mat>
class ParticleTracker{
public:

	ParticleTracker(std::string filename, const std::vector<Particle<Vec,Mat>*>& parts):
		parts_(parts){
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


	virtual void output() = 0;
protected:
	std::ostream* out_;
	const std::vector<Particle<Vec,Mat>* >& parts_;
};

}//output namespace
}//treecode namespace

#endif /* PARTICLETRACKER_H_ */
