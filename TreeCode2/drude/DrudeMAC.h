/*
 * DrudeMAC.h
 *
 *  Created on: 10 Mar 2012
 *      Author: stefans
 */

#ifndef DRUDEMAC_H_
#define DRUDEMAC_H_

#include "../macs/AcceptanceCriterion.h"
#include "../bounds/BoundaryConditions.h"
#include "BigParticle.h"

namespace treecode{
template <class Vec, class Mat>
class DrudeMAC : public AcceptanceCriterion<Vec, Mat>{
public:
	DrudeMAC(const BoundaryConditions<Vec,Mat>& bounds, double dt):bounds_(bounds),dt_(dt){}

	typename AcceptanceCriterion<Vec,Mat>::result accept(const Particle<Vec,Mat>& p, const Node<Vec, Mat>& n) const {
		if(n.getStatus() != Node<Vec,Mat>::LEAF)
			return AcceptanceCriterion<Vec,Mat>::CONTINUE;

		Particle<Vec,Mat>* node_part = n.getParticles().front();
		BigParticle<Vec,Mat>* bp = dynamic_cast<BigParticle<Vec,Mat>*>(node_part);
		if(bp == NULL)
			return AcceptanceCriterion<Vec,Mat>::CONTINUE;

		//If the next step will take us into a particle, then interact with that node.
		if(DrudeMAC<Vec,Mat>::willIntersect(p, *bp, bounds_)){
			double dist = DrudeMAC<Vec,Mat>::distanceToIntersection(p, *bp, bounds_);
			double will_travel = p.getVelocity().norm() * dt_;

			if(dist > 0 && dist < will_travel)
				return AcceptanceCriterion<Vec,Mat>::ACCEPT;
			else
				return AcceptanceCriterion<Vec,Mat>::REJECT;
		}else
			return AcceptanceCriterion<Vec,Mat>::REJECT;
	}

	void setTimestep(double dt){dt_ = dt;}

	static bool willIntersect(const Particle<Vec,Mat>& p, const BigParticle<Vec,Mat>& bp, const BoundaryConditions<Vec,Mat>& bc){
		Vec l = p.getVelocity() / p.getVelocity().norm();
		Vec c = bc.getDisplacementVector(bp.getPosition(), p.getPosition());
		return (l.dot(c) * l.dot(c) - c.squaredNorm() + bp.getRadius()*bp.getRadius()) >= 0;
	}

	static double distanceToIntersection(const Particle<Vec,Mat>& p, const BigParticle<Vec,Mat>& bp, const BoundaryConditions<Vec,Mat>& bc){
		Vec l = p.getVelocity() / p.getVelocity().norm();
		Vec c = bc.getDisplacementVector(bp.getPosition(), p.getPosition());
		double determinant = l.dot(c) * l.dot(c) - c.squaredNorm() + bp.getRadius()*bp.getRadius();
		if(determinant < 0){
			std::cout << " ";
		}
		return l.dot(c) - sqrt(determinant);
	}

private:
	const BoundaryConditions<Vec,Mat>& bounds_;
	double dt_;
};
}

#endif /* DRUDEMAC_H_ */
