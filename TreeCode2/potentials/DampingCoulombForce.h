/*
 * CoulombForceEField.h
 *
 *  Created on: 20 Feb 2012
 *      Author: stefans
 */

#ifndef DAMPING_COULOMBFORCE_
#define DAMPING_COULOMBFORCE_

#include "CoulombForce.h"
#include "../bounds/BoundaryConditions.h"

namespace treecode {
namespace potentials {

template <int D>
class DampingCoulombForce : public CoulombForceThreeD<D>{
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;

public:
	DampingCoulombForce(double force_softening, const BoundaryConditions<D>& bc, double damping):
		CoulombForceThreeD<D>(force_softening, bc), damping_(damping){}

	virtual Vec getForce(const Particle<D>& part, const Node<D>& node, Precision precision) const{
		Vec force = CoulombForceThreeD<D>::getForce(part, node, precision);
		force -= damping_ * part.getVelocity();
		return force;
	}

	virtual ~DampingCoulombForce(){}
private:
	double damping_;
};

} /* namespace output */
} /* namespace treecode */
#endif /* DAMPING_COULOMBFORCE_ */
