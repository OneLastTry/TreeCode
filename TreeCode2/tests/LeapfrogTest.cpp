/*
 * TreeTest.cpp
 *
 *  Created on: 6 Dec 2011
 *      Author: stefans
 */


#include <boost/test/unit_test.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <cmath>
#include <Eigen/Dense>

#include <iostream>

#include "../Tree.h"
#include "../Node.h"
#include "../Particle.h"
#include "../Configuration.h"
#include "../bounds/OpenBoundary.h"
#include "../potentials/CoulombForce.h"
#include "../pushers/LeapfrogPusher.h"

using namespace std;
using namespace treecode;
using namespace Eigen;
using namespace treecode::pusher;
using namespace treecode::potentials;

#define TOLERANCE 1E-5
#define MAX_TIME 10

#ifdef __CDT_PARSER__
#define BOOST_CHECK_CLOSE(a,b,c)
#endif

BOOST_AUTO_TEST_SUITE(PusherTests)

BOOST_AUTO_TEST_CASE(CoulombPotentialTest2D){
	boost::random::mt19937 rng;
	boost::random::uniform_real_distribution<double> rand_dist(2, 10);

	Configuration c(2, 0.0, 0.01, 0.01, 1.0/3);

	//Create a particle at the orign
	Particle p1(c, +1, 1, Vector2d(0, 0), Vector2d(0,0));
	vector<Particle*> p;
	p.push_back(&p1);

	//Manually form a node around p1
	Node n1(c, Vector2d(-1,-1), 2);
	n1.setParticles(p);
	n1.setStatus(Node::tree_status::LEAF);
	n1.calculateMonopoleMoment();
	n1.calculateDipoleMoment();
	n1.calculateQuadrupoleMoment();

	//Init a force
	OpenBoundary ob(c);
	ob.init(p);
	CoulombForceThreeD cf(c, ob);

	//Create a randomly placed test particle, and test that the directly
	//calculated force is equal.
	Particle p2(c, +1, 1, Vector2d(rand_dist(rng), rand_dist(rng)), Vector2d::Zero());

	Vector2d node_force = cf.getForce(p2, n1, Potential::Precision::quadrupole);
	Vector2d disp_vec = p2.getPosition() - p1.getPosition();	//r
	Vector2d anal_force = p1.getCharge() * p2.getCharge() * disp_vec
			/ disp_vec.squaredNorm() / disp_vec.norm();	//q1 q1 r / |r|^3

	BOOST_CHECK_CLOSE(anal_force[0], node_force[0], TOLERANCE);
	BOOST_CHECK_CLOSE(anal_force[1], node_force[1], TOLERANCE);

	//Do the same with the potential
	double node_potential = cf.getPotential(p2, n1, Potential::Precision::quadrupole);
	double anal_potential = p1.getCharge() * p2.getCharge() / disp_vec.norm();
	BOOST_CHECK_CLOSE(anal_potential, node_potential, TOLERANCE);

}

/**
 * Check that the leapfrog is time-symmetric.
 */
BOOST_AUTO_TEST_CASE(LFPusher){
	boost::random::mt19937 rng;
	boost::random::uniform_real_distribution<double> uniform_dist(0, 20);
	boost::random::bernoulli_distribution<double> boolean_dist(0.5);

	Configuration c(2, 0.0, 0.01, 0.01, 1.0/3);

	vector<VectorXd> original_positions;

	//Setup some test particles
	vector<Particle*> leapfrog_parts;
	leapfrog_parts.push_back(new Particle(c, +1, 10000, VectorXd::Zero(c.getNumDimensions()), VectorXd::Zero(c.getNumDimensions())));
	VectorXd pos(c.getNumDimensions());	pos << 1,0;
	VectorXd vel(c.getNumDimensions());	vel << 0,1;
	leapfrog_parts.push_back(new Particle(c, -1, 1, pos, vel));


	OpenBoundary ob(c);
	ob.init(leapfrog_parts);
	CoulombForceThreeD cf(c, ob);

	Tree tr(c, ob, leapfrog_parts);

	LeapfrogPusher lp(c, ob, cf);
	lp.init(leapfrog_parts, tr, Potential::Precision::quadrupole);

	int num_timesteps = MAX_TIME / c.getTimestep();
	for(int i = 0; i < num_timesteps; i++){
		tr.rebuild();
		//Push all particles
		lp.push_particles(leapfrog_parts, tr, ob, Potential::Precision::quadrupole);
		ob.timestepOver();
	}
	//Reverse time.
	c.setTimestep(-1 * c.getTimestep());
	for(int i = 0; i < num_timesteps; i++){
		tr.rebuild();
		//Push all particles
		lp.push_particles(leapfrog_parts, tr, ob, Potential::Precision::quadrupole);
		ob.timestepOver();
	}

	for(unsigned int i=0;i<original_positions.size();i++){
		for(unsigned int j=0;j<c.getNumDimensions();j++){
			//Check to 0.1%
			BOOST_CHECK_CLOSE(leapfrog_parts[i]->getPosition()[j], original_positions[i][j], 0.1);
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
