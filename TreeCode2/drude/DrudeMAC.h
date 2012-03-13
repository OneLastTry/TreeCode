/*
 * DrudeMAC.h
 *
 *  Created on: 10 Mar 2012
 *      Author: stefans
 */

#ifndef DRUDEMAC_H_
#define DRUDEMAC_H_

#include "../macs/AcceptanceCriterion.h"
#include "BigParticle.h"

namespace treecode{
template <class Vec, class Mat>
class DrudeMAC : public AcceptanceCriterion<Vec, Mat>{
public:
	DrudeMAC(double dt):dt_(dt){}

	bool accept(const Particle<Vec>& p, const Node<Vec, Mat>& n) const {
		if(n.getStatus() != Node<Vec,Mat>::LEAF)
			return false;
		Particle<Vec>* node_part = n.getParticles().front();
		BigParticle<Vec>* bp = dynamic_cast<BigParticle<Vec>*>(node_part);
		if(bp == NULL)
			return false;

		//If the next step will take us into a particle, then interact with that node.
		Vec next_pos = p.getPosition() + p.getVelocity()*dt_;
		Vec disp_vec = next_pos - n.getCentreOfCharge();

		return (disp_vec.norm() < bp->getRadius());
	}
private:
	double dt_;
};
}

#endif /* DRUDEMAC_H_ */
