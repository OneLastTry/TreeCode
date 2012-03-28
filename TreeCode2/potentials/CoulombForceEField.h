/*
 * CoulombForceEField.h
 *
 *  Created on: 20 Feb 2012
 *      Author: stefans
 */

#ifndef COULOMBFORCEEFIELD_H_
#define COULOMBFORCEEFIELD_H_

#include "CoulombForce.h"
#include "../bounds/BoundaryConditions.h"

namespace treecode {
namespace potentials {

template <int D>
class CoulombForceEField : public CoulombForceThreeD<D>{
	typedef Eigen::Matrix<double, D, D> Mat;
	typedef Eigen::Matrix<double, D, 1> Vec;

public:
	CoulombForceEField(double force_softening, const BoundaryConditions<D>& bc, Vec e_field):
		CoulombForceThreeD<D>(force_softening, bc), e_field_(e_field){}

	virtual Vec getForce(const Particle<D>& part, const Node<D>& node, Precision precision) const{
		Vec force = CoulombForceThreeD<D>::getForce(part, node, precision);
		force += part.getCharge() * e_field_;
		return force;
	}

	virtual ~CoulombForceEField(){}
private:
	Vec e_field_;
};

} /* namespace output */
} /* namespace treecode */
#endif /* COULOMBFORCEEFIELD_H_ */
