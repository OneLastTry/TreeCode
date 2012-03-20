/*
 * TreeTest.cpp
 *
 *  Created on: 6 Dec 2011
 *      Author: stefans
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <cmath>

#include <Eigen/Dense>
#include <Particle.h>
#include <distributions/UniformDistribution.h>
#include <distributions/ConstDistribution.h>
#include <distributions/ConstantChargeDistribution.h>


using namespace std;
using namespace treecode;
using namespace Eigen;
using namespace treecode::distribution;

#define TOLERANCE 1E-5
#define FORCE_TOLERANCE 0.05

template<typename V, typename M>
struct F {
    F(){
    	using boost::random::mt19937;
    	mt19937 rng;

    	double length, temperature;
    	unsigned int num_particles;

    	UniformDistribution3d 			position_dist(Vec::Zero(), Vec::Zero().array() + length);
    	ConstDistribution3d				velocity_dist(Vec::Zero());
    	ConstantChargeDistribution3d	electron_charges(-1);
    	ConstantChargeDistribution3d	ion_charges(1);

    	int id = 0;
    	vector<Particle3d*> ions = Particle3d::generateParticles<mt19937>(num_particles, 1837, rng,
    			position_dist, velocity_dist, ion_charges, id);
    	vector<Particle3d*> electrons = Particle3d::generateParticles<mt19937>(num_particles, 1837, rng,
    	    			position_dist, velocity_dist, ion_charges, id);
    	parts.insert(parts.end(), ions.begin(), ions.end());
    	parts.insert(parts.end(), electrons.begin(), electrons.end());
    	OpenBoundary<V,M> bounds();
    	tree = new Tree()
    }

    ~F(){
    	delete tree;
    }

    Tree<V,M> *tree;
    std::vector<Particle<V,M>* > parts;
};

template<typename V, typename M>
V analytic_dipole_moment(const vector<Particle<V,M>*>& parts, const V& origin){
	typedef Particle<V,M> Particle;
	V dp_moments = V::Zero();
	BOOST_FOREACH(Particle *p, parts){
		dp_moments += p->getCharge() * (p->getPosition() - origin);
	}
	return dp_moments;
}

template<typename V, typename M>
M analytic_quadrupole_moment(const vector<Particle<V,M>*> parts, const V& origin){
	typedef Particle<V,M> Particle;
	M qm = M::Zero();
	BOOST_FOREACH(Particle *p, parts){
		V disp_vec = p->getPosition() - origin;
		qm += (disp_vec * disp_vec.transpose()) * p->getCharge();
	}
	return qm;
}

template<typename V, typename M>
V analytic_centre_of_charge(const vector<Particle<V,M>*> parts, int& charge, int& abs_charge){
	charge = 0;
	abs_charge = 0;
	V coc = V::Zero();
	typedef Particle<V,M> Particle;
	BOOST_FOREACH(Particle* p, parts){
		coc += abs(p->getCharge()) * p->getPosition();
		charge += p->getCharge();
		abs_charge += abs(p->getCharge());
	}
	coc /= abs_charge;
	return coc;
}
template<typename V, typename M>
V analytic_force(const Particle<V,M>& test_part, const vector<Particle<V,M>*> parts){
	V force = V::Zero();
	typedef Particle<V,M> Particle;
	BOOST_FOREACH(Particle* p, parts){
		V r = (test_part.getPosition() - p->getPosition());
		force += p->getCharge() * r / r.squaredNorm() / r.norm();
	}
	return force;
}
template<typename V, typename M>
double analytic_potential(const Particle<V,M>& test_part, const vector<Particle<V,M>*> parts){
	double potential = 0;
	typedef Particle<V,M> Particle;
	BOOST_FOREACH(Particle* p, parts){
		V r = (test_part.getPosition() - p->getPosition());
		potential += p->getCharge() / r.norm();
	}
	return potential;
}

BOOST_AUTO_TEST_SUITE(Moments);

BOOST_AUTO_TEST_CASE(TwoDMoments){

}

BOOST_AUTO_TEST_SUITE_END();
