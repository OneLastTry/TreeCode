/*
 * analytic.cpp
 *
 *  Created on: 30 Mar 2012
 *      Author: stefans
 */

#ifndef _ANALYTIC_H_
#define _ANALYTIC_H_

#include <boost/foreach.hpp>
#include <Eigen/Dense>
#include <vector>
#include <cmath>

#include <Tree.h>
#include <Particle.h>

namespace treecode{
namespace tests{

template <int D>
class AnalyticResults {
	typedef Eigen::Matrix<double, D, D> M;
	typedef Eigen::Matrix<double, D, 1> V;
public:

V analytic_dipole_moment(const std::vector<Particle<D>*>& parts, const V& origin){
	typedef Particle<D> Particle;
	V dp_moments = V::Zero();
	BOOST_FOREACH(Particle *p, parts){
		dp_moments += p->getCharge() * (p->getPosition() - origin);
	}
	return dp_moments;
}

M analytic_quadrupole_moment(const std::vector<Particle<D>*> parts, const V& origin){
	typedef Particle<D> Particle;
	M qm = M::Zero();
	BOOST_FOREACH(Particle *p, parts){
		V disp_vec = p->getPosition() - origin;
		qm += (disp_vec * disp_vec.transpose()) * p->getCharge();
	}
	return qm;
}

V analytic_centre_of_charge(const std::vector<Particle<D>*> parts, int& charge, int& abs_charge){
	charge = 0;
	abs_charge = 0;
	V coc = V::Zero();
	typedef Particle<D> Particle;
	BOOST_FOREACH(Particle* p, parts){
		coc += abs(p->getCharge()) * p->getPosition();
		charge += p->getCharge();
		abs_charge += abs(p->getCharge());
	}
	coc /= abs_charge;
	return coc;
}

V analytic_force(const Particle<D>& test_part, const std::vector<Particle<D>*> parts){
	V force = V::Zero();
	typedef Particle<D> Particle;
	BOOST_FOREACH(Particle* p, parts){
		if(p == &test_part)
			continue;
		V r = (test_part.getPosition() - p->getPosition());
		force += p->getCharge() * r / r.squaredNorm() / r.norm();
	}
	return force;
}

double analytic_potential(const Particle<D>& test_part, const std::vector<Particle<D>*> parts){
	double potential = 0;
	typedef Particle<D> Particle;
	BOOST_FOREACH(Particle* p, parts){
		if(p == &test_part)
			continue;
		V r = (test_part.getPosition() - p->getPosition());
		potential += p->getCharge() / r.norm();
	}
	return potential;
}

};

}
}

#endif /* _ANALYTIC_H_ */
