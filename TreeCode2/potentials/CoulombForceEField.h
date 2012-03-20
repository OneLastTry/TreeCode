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

template <class Vec, class Mat>
class CoulombForceEField : public CoulombForceThreeD<Vec, Mat>{
public:
	CoulombForceEField(double force_softening, const BoundaryConditions<Vec,Mat>& bc, Vec e_field):
		CoulombForceThreeD<Vec, Mat>(force_softening, bc), e_field_(e_field){}

	virtual Vec getForce(const Particle<Vec,Mat>& part, const Node<Vec,Mat>& node, Precision precision) const{
		Vec force = CoulombForceThreeD<Vec, Mat>::getForce(part, node, precision);
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
