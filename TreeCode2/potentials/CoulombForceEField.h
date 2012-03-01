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
#include "../Configuration.h"

namespace treecode {
namespace potentials {

template <class Vec, class Mat>
class CoulombForceEField : public CoulombForceThreeD<Vec, Mat>{
public:
	CoulombForceEField(const Configuration<Vec>& conf, const BoundaryConditions<Vec>& bc, Vec e_field):
		CoulombForceThreeD<Vec, Mat>(conf, bc), e_field_(e_field){}

	virtual Vec getForce(const Particle<Vec>& part, const Node<Vec,Mat>& node, Precision precision) const{
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
