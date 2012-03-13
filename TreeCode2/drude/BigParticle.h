/*
 * DrudeParticle.h
 *
 *  Created on: 10 Mar 2012
 *      Author: stefans
 */

#ifndef BIGPARTICLE_H_
#define BIGPARTICLE_H_

#include "../Particle.h"

namespace treecode{
template <class Vec>
class BigParticle : public Particle<Vec>{
public:
	BigParticle(int q, int m, const Vec& pos, const Vec& vel, unsigned int id, double radius):
		Particle<Vec>(q,m,pos,vel,id), radius_(radius){}

	double getRadius(){return radius_;}

private:
	double radius_;
};
}

#endif /* BIGPARTICLE_H_ */
