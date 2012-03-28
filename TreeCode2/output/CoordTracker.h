/*
 * ParticleTracker.h
 *
 *  Created on: 3 Mar 2012
 *      Author: stefans
 */

#ifndef COORDTRACKER_H_
#define COORDTRACKER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/foreach.hpp>
#include "../Particle.h"
#include "ParticleTracker.h"

namespace treecode{
namespace output{

template <class Vec, class Mat>
class CoordTracker : public ParticleTracker<Vec,Mat>{
public:
	enum token {
		POSITION,
		VELOCITY
	};

	CoordTracker(std::string filename, const std::vector<Particle<Vec,Mat>*>& parts, token record):
		ParticleTracker<Vec,Mat>(filename, parts), record_(record){	}

	virtual ~CoordTracker(){}


	virtual void output(){
		typedef ParticleTracker<Vec,Mat> parent;
		typedef Particle<Vec,Mat> part_t;
		BOOST_FOREACH(part_t* p, parent::parts_){
			if(record_ == POSITION){
				for(unsigned int i = 0; i < p->getPosition().rows();i++)
					(*(parent::out_)) << p->getPosition()[i] << "\t";
			}else if(record_ == VELOCITY){
				for(unsigned int i = 0; i < p->getVelocity().rows();i++)
					(*(parent::out_)) << p->getVelocity()[i] << "\t";
			}
		}
		(*(parent::out_)) << std::endl;
	}
protected:
	token record_;
};

}//output namespace
}//treecode namespace

#endif /* COORDTRACKER_H_ */
