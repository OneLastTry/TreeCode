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
template <int D>
class BarnesHutMAC : public AcceptanceCriterion<D>{
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;
public:
	BarnesHutMAC(double theta, const BoundaryConditions<D>& bounds):theta_(theta), bounds_(bounds){}
	typename AcceptanceCriterion<D>::result accept(const Particle<D>& p, const Node<D>& n) const{
        double d_squared = bounds_.getDisplacementVector(p.getPosition(), n.getCentreOfCharge()).squaredNorm();
        double mac = n.getSize() * n.getSize() / d_squared;
        //Either add to ilist or recurse into daughters.
        if(mac < theta_*theta_ || n.getStatus() == Node<D>::LEAF)
        	return AcceptanceCriterion<D>::ACCEPT;
        else
        	return AcceptanceCriterion<D>::CONTINUE;
	}

	double getTheta() const {return theta_;}

	void setTheta(double theta){
		theta_ = theta;
	}

private:
	double theta_;
	const BoundaryConditions<D>& bounds_;
};
}

#endif /* BARNESHUTMAC_H_ */
