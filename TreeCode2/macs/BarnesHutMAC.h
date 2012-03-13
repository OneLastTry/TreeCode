/*
 * BarnesHutMAC.h
 *
 *  Created on: 10 Mar 2012
 *      Author: stefans
 */

#ifndef BARNESHUTMAC_H_
#define BARNESHUTMAC_H_

#include "AcceptanceCriterion.h"
#include "../bounds/BoundaryConditions.h"

namespace treecode{
template <class Vec, class Mat>
class BarnesHutMAC : public AcceptanceCriterion<Vec,Mat>{
public:
	BarnesHutMAC(double theta, const BoundaryConditions<Vec>& bounds):theta_(theta), bounds_(bounds){}
	bool accept(const Particle<Vec>& p, const Node<Vec,Mat>& n) const{
        double d_squared = bounds_.getDisplacementVector(p.getPosition(), n.getCentreOfCharge()).squaredNorm();
        double mac = n.getSize() * n.getSize() / d_squared;
        //Either add to ilist or recurse into daughters.
        return (mac < theta_*theta_ || n.getStatus() == Node<Vec,Mat>::LEAF);
	}

	double getTheta() const {return theta_;}

	void setTheta(double theta){
		theta_ = theta;
	}

private:
	double theta_;
	const BoundaryConditions<Vec>& bounds_;
};
}

#endif /* BARNESHUTMAC_H_ */
