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
template <int D>
class DrudeMAC : public AcceptanceCriterion<D>{
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;

public:
	DrudeMAC(const BoundaryConditions<D>& bounds, double dt):bounds_(bounds),dt_(dt){}

	typename AcceptanceCriterion<D>::result accept(const Particle<D>& p, const Node<D>& n) const {
		if(n.getStatus() != Node<D>::LEAF)
			return AcceptanceCriterion<D>::CONTINUE;

		Particle<D>* node_part = n.getParticles().front();
		BigParticle<D>* bp = dynamic_cast<BigParticle<D>*>(node_part);
		if(bp == NULL)
			return AcceptanceCriterion<D>::CONTINUE;

		//If the next step will take us into a particle, then interact with that node.
		if(DrudeMAC<D>::willIntersect(p, *bp, bounds_)){
			double dist = DrudeMAC<D>::distanceToIntersection(p, *bp, bounds_);
			double will_travel = p.getVelocity().norm() * dt_;

			if(dist > 0 && dist < will_travel)
				return AcceptanceCriterion<D>::ACCEPT;
			else
				return AcceptanceCriterion<D>::REJECT;
		}else
			return AcceptanceCriterion<D>::REJECT;
	}

	void setTimestep(double dt){dt_ = dt;}

	static bool willIntersect(const Particle<D>& p, const BigParticle<D>& bp, const BoundaryConditions<D>& bc){
		Vec l = p.getVelocity() / p.getVelocity().norm();
		Vec c = bc.getDisplacementVector(bp.getPosition(), p.getPosition());
		return (l.dot(c) * l.dot(c) - c.squaredNorm() + bp.getRadius()*bp.getRadius()) >= 0;
	}

	static double distanceToIntersection(const Particle<D>& p, const BigParticle<D>& bp, const BoundaryConditions<D>& bc){
		Vec l = p.getVelocity() / p.getVelocity().norm();
		Vec c = bc.getDisplacementVector(bp.getPosition(), p.getPosition());
		double determinant = l.dot(c) * l.dot(c) - c.squaredNorm() + bp.getRadius()*bp.getRadius();
		if(determinant < 0){
			std::cout << " ";
		}
		return l.dot(c) - sqrt(determinant);
	}

private:
	const BoundaryConditions<D>& bounds_;
	double dt_;
};
}

#endif /* DRUDEMAC_H_ */
