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
template <class Vec, class Mat>
class BigParticle : public Particle<Vec,Mat>{
public:
	BigParticle(int q, int m, const Vec& pos, const Vec& vel, unsigned int id, double radius):
		Particle<Vec,Mat>(q,m,pos,vel,id), radius_(radius){}

	double getRadius()const {return radius_;}

private:
	double radius_;
};
}

#endif /* BIGPARTICLE_H_ */
