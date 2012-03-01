/*
 * TreeTest.cpp
 *
 *  Created on: 6 Dec 2011
 *      Author: stefans
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <cmath>

#include <Eigen/Dense>
#include "../Tree.h"
#include "../Node.h"
#include "../Particle.h"
#include "../Configuration.h"
#include "../bounds/OpenBoundary.h"
#include "../potentials/CoulombForce.h"

using namespace std;
using namespace treecode;
using namespace Eigen;

#define TOLERANCE 1E-5
#define FORCE_TOLERANCE 0.05

#ifdef __CDT_PARSER__
#define BOOST_CHECK_CLOSE(a,b,c)
#endif

MatrixXd analytic_quadrupole_moment(const Configuration& c, const vector<Particle*> parts, const VectorXd& origin){
	MatrixXd qm = MatrixXd::Zero(c.getNumDimensions(), c.getNumDimensions());
	for(Particle* p : parts){
		VectorXd disp_vec = p->getPosition() - origin;
		qm += (disp_vec * disp_vec.transpose()) * p->getCharge();
	}
	return qm;
}

VectorXd analytic_dipole_moment(const Configuration& c, const vector<Particle*>& parts, const VectorXd& origin){
	VectorXd dp_moments = VectorXd::Zero(c.getNumDimensions());
	for(Particle* p : parts)
		dp_moments += p->getCharge() * (p->getPosition() - origin);
	return dp_moments;
}

VectorXd analytic_centre_of_charge(const Configuration& c, const vector<Particle*> parts, int& charge, int& abs_charge){
	charge = 0;
	abs_charge = 0;
	VectorXd coc = VectorXd::Zero(c.getNumDimensions());
	for(Particle* p : parts){
		coc += abs(p->getCharge()) * p->getPosition();
		charge += p->getCharge();
		abs_charge += abs(p->getCharge());
	}
	coc /= abs_charge;
	return coc;
}

VectorXd analytic_force(const Configuration& c, const Particle& test_part, const vector<Particle*> parts){
	VectorXd force = VectorXd::Zero(c.getNumDimensions());
	for(Particle* p : parts){
		VectorXd r = (test_part.getPosition() - p->getPosition());
		force += p->getCharge() * r / r.squaredNorm() / r.norm();
	}
	return force;
}

double analytic_potential(const Configuration& c, const Particle& test_part, const vector<Particle*> parts){
	double potential = 0;
	for(Particle* p : parts){
		VectorXd r = (test_part.getPosition() - p->getPosition());
		potential += p->getCharge() / r.norm();
	}
	return potential;
}

BOOST_AUTO_TEST_SUITE(Moments)

BOOST_AUTO_TEST_CASE(TwoDMoments){
	boost::random::mt19937 rng;
	boost::random::uniform_01<double> uniform_dist;
	boost::random::bernoulli_distribution<bool> boolean_dist(0.5);

	Configuration c(2, 0.0, 0.01, 0.01, 1.0/3);

	//Setup some test particles
	vector<Particle*> parts;
	for (int i = 0; i < 2; i++) {
		VectorXd pos(c.getNumDimensions());
		for(unsigned int j=0;j<c.getNumDimensions();j++)
			pos[j] = uniform_dist(rng);
		int charge = (uniform_dist(rng) > 0.5) ? -1 : +1;
		parts.push_back( new Particle(c, charge, 1, pos, VectorXd::Zero(3)) );
	}

	OpenBoundary ob(c);
	ob.init(parts);

	Tree tr(c, ob, parts);
	tr.rebuild();

	int observed_charge, observed_abs_charge;
	VectorXd observed_coc = analytic_centre_of_charge(c, parts, observed_charge, observed_abs_charge);
	VectorXd observed_dipole_moments = analytic_dipole_moment(c, parts, observed_coc);
	MatrixXd observed_quadrupole_moments = analytic_quadrupole_moment(c, parts, observed_coc);

	int tree_charge = tr.getTotalCharge();
	int tree_abs_charge = tr.getTotalAbsCharge();
	VectorXd tree_coc = tr.getTotalCentreOfCharge();
	VectorXd tree_dp = tr.getTotalDipoleMoments();
	MatrixXd tree_qp = tr.getTotalQuadrupoleMoments();

	BOOST_CHECK_EQUAL(observed_charge, tree_charge);
	BOOST_CHECK_EQUAL(observed_abs_charge, tree_abs_charge);

	//Check each component of dipole moment
	for(unsigned int i=0;i<c.getNumDimensions();i++){
		BOOST_CHECK_CLOSE(observed_coc[i], tree_coc[i], TOLERANCE);
		BOOST_CHECK_CLOSE(observed_dipole_moments[i], tree_dp[i], TOLERANCE);
		for (unsigned int j = 0; j < c.getNumDimensions(); j++) {
			BOOST_CHECK_CLOSE(observed_quadrupole_moments(i,j), tree_qp(i,j), TOLERANCE);
		}
	}

	VectorXd test_pos = VectorXd::Zero(c.getNumDimensions());
	test_pos[0] = 100;
	Particle test_particle(c, 1, 1, test_pos, VectorXd::Zero(c.getNumDimensions()));
	potentials::CoulombForceThreeD cf(c, ob);

	VectorXd direct_force = analytic_force(c, test_particle, parts);
	VectorXd tree_force = cf.getForce(test_particle, tr.getRoot(), potentials::Potential::Precision::quadrupole);

	double direct_potential = analytic_potential(c, test_particle, parts);
	double tree_potential = cf.getPotential(test_particle, tr.getRoot(), potentials::Potential::Precision::quadrupole);

	BOOST_CHECK_CLOSE(direct_potential, tree_potential, FORCE_TOLERANCE);
	for(unsigned int i=0;i<c.getNumDimensions();i++){
		BOOST_CHECK_CLOSE(direct_force[i], tree_force[i], FORCE_TOLERANCE);
	}

}

BOOST_AUTO_TEST_SUITE_END()
