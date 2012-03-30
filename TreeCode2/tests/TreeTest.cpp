/*
 * TreeTest.cpp
 *
 *  Created on: 6 Dec 2011
 *      Author: stefans
 */


#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

#include <Particle.h>
#include <distributions/UniformDistribution.h>
#include <distributions/ConstDistribution.h>
#include <distributions/ConstantChargeDistribution.h>
#include <Tree.h>
#include <Node.h>
#include <bounds/OpenBoundary.h>
#include <potentials/CoulombForce.h>
#include <potentials/Potential.h>
#include <macs/BarnesHutMAC.h>

#include "analytic.h"
#include "custom_asserts.h"

#define TOLERANCE 1E-5
#define FORCE_TOLERANCE 0.05

/**
 * @brief Fixture for tree tests.
 */
template<int D>
struct TreeFixture {
	typedef Eigen::Matrix<double, D, 1> V;

	/**
	 * Create a list of 1000 particles, and create a tree.
	 */
    TreeFixture(){
    	using boost::random::mt19937;
    	using treecode::distribution::UniformDistribution;
    	using treecode::distribution::ConstDistribution;
    	using treecode::distribution::ConstantChargeDistribution;
    	using treecode::Particle;

    	mt19937 rng;

    	double length = 1.0;
    	unsigned int num_particles = 1000;

    	//Create distributions
    	UniformDistribution<mt19937> 	position_dist(rng, V::Zero(), V::Zero().array() + length);
    	ConstDistribution				velocity_dist(V::Zero());
    	ConstantChargeDistribution		electron_charges(-1);
    	ConstantChargeDistribution		ion_charges(1);

    	//Generate some particles
    	std::vector<Particle<D>*> ions = Particle<D>::template generateParticles<mt19937>(num_particles, 1837, rng,
    			position_dist, velocity_dist, ion_charges);
    	std::vector<Particle<D>*> electrons = Particle<D>::template generateParticles<mt19937>(num_particles, 1, rng,
    	    			position_dist, velocity_dist, electron_charges);
    	parts.insert(parts.end(), ions.begin(), ions.end());
    	parts.insert(parts.end(), electrons.begin(), electrons.end());
    	bounds.init(parts);
    	//Create a tree
    	tree = new treecode::Tree<D>(bounds, parts);
    	tree->rebuild();
    }

    //Get force on a single particle from tree
    V force_on_particle(
    		treecode::Particle<D>* p,
    		const treecode::potentials::Potential<D>& potential,
    		treecode::potentials::Precision precision,
    		const treecode::AcceptanceCriterion<D>& mac
    		) {
    	typedef treecode::Node<D> Node;
    	typedef std::vector<Node*> interaction_list;
    	V force = V::Zero();

    	//Build interaction list, and add force
    	interaction_list ilist;
    	tree->getInteractionList(*p, ilist, mac);
    	for (typename interaction_list::iterator it = ilist.begin(); it < ilist.end(); it++) {
    		Node* n = *it;
    		force += potential.getForce(*p, *n, precision);
    	}
    	return force;
    }

    //Get potential at particle from tree
    double potential_at_particle(
    		treecode::Particle<D>* p,
    		const treecode::potentials::Potential<D>& potential,
    		treecode::potentials::Precision precision,
    		const treecode::AcceptanceCriterion<D>& mac
    		) {
    	typedef treecode::Node<D> Node;
    	typedef std::vector<Node*> interaction_list;
    	double pot = 0;

    	//Build interaction list and get potential
    	interaction_list ilist;
    	tree->getInteractionList(*p, ilist, mac);
    	for (typename interaction_list::iterator it = ilist.begin(); it < ilist.end(); it++) {
    		Node* n = *it;
    		pot += potential.getPotential(*p, *n, precision);
    	}
    	return pot;
    }

    ~TreeFixture(){
    	treecode::Particle<D>::deleteParticles(parts);
    	delete tree;
    }

    treecode::OpenBoundary<D> bounds;
    treecode::Tree<D> *tree;
    std::vector<treecode::Particle<D>* > parts;

    //Analytic results
    treecode::tests::AnalyticResults<D> anal;
};



typedef Eigen::Vector2d V2d;
typedef Eigen::Matrix2d M2d;
BOOST_FIXTURE_TEST_SUITE(Moments2d, TreeFixture<2>)

//Compare analytic coc, dipole and quadrupole moments.
BOOST_AUTO_TEST_CASE(TwoDMoments){
	const treecode::Node<2>& root = tree->getRoot();
	int total_charge, modulus_charge;
	V2d anal_coc  = anal.analytic_centre_of_charge(parts, total_charge, modulus_charge);
	V2d anal_dip  = anal.analytic_dipole_moment(parts, root.getCentreOfCharge());
	M2d anal_quad = anal.analytic_quadrupole_moment(parts, root.getCentreOfCharge());

	EIGEN_REQUIRE_CLOSE(anal_coc, root.getCentreOfCharge(), TOLERANCE);
	EIGEN_REQUIRE_CLOSE(anal_dip, root.getDipoleMoments(), TOLERANCE);
	EIGEN_REQUIRE_CLOSE(anal_quad, root.getQuadrupoleMoments(), TOLERANCE);
}

//Compare analytic forces, and check forces from tree
BOOST_AUTO_TEST_CASE(TwoDForces){
	treecode::potentials::CoulombForceThreeD<2> pot(0.0, bounds);
	treecode::BarnesHutMAC<2> mac(0.0, bounds);
	treecode::Particle<2>* part = parts.front();

	V2d tree_force = force_on_particle(part, pot, treecode::potentials::quadrupole, mac);
	V2d anal_force = anal.analytic_force(*part, parts);
	double tree_pot = potential_at_particle(part, pot, treecode::potentials::quadrupole, mac);
	double anal_pot = anal.analytic_potential(*part, parts);

	//Should be identical with theta = 0
	EIGEN_REQUIRE_CLOSE(anal_force, tree_force, TOLERANCE);
	BOOST_REQUIRE_CLOSE(tree_pot, anal_pot, TOLERANCE);

	//Now try with theta = 0.5
	mac.setTheta(0.5);
	V2d mp_force = force_on_particle(part, pot, treecode::potentials::monopole, mac);
	double mp_pot = potential_at_particle(part, pot, treecode::potentials::monopole, mac);

	V2d dp_force = force_on_particle(part, pot, treecode::potentials::dipole, mac);
	double dp_pot = potential_at_particle(part, pot, treecode::potentials::dipole, mac);

	V2d qp_force = force_on_particle(part, pot, treecode::potentials::quadrupole, mac);
	double qp_pot = potential_at_particle(part, pot, treecode::potentials::quadrupole, mac);

	BOOST_TEST_MESSAGE("Is the following error acceptable? (theta = 0.5)");
	BOOST_TEST_MESSAGE("Analytic force:   (" << anal_force[0] << ", " << anal_force[1] << ")");
	BOOST_TEST_MESSAGE("Monopole force:   (" << mp_force[0] << ", " << mp_force[1] << ")");
	BOOST_TEST_MESSAGE("Dipole force:     (" << dp_force[0] << ", " << dp_force[1] << ")");
	BOOST_TEST_MESSAGE("Quadrupole force: (" << qp_force[0] << ", " << qp_force[1] << ")");
	BOOST_TEST_MESSAGE("Analytic potential:   " << anal_pot);
	BOOST_TEST_MESSAGE("Monopole potential:   " << mp_pot);
	BOOST_TEST_MESSAGE("Dipole potential:     " << dp_pot);
	BOOST_TEST_MESSAGE("Quadrupole potential: " << qp_pot);
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_FIXTURE_TEST_SUITE(Moments3d, TreeFixture<3>)

typedef Eigen::Vector3d V3d;
typedef Eigen::Matrix3d M3d;
//Compare analytic coc, dipole and quadrupole moments.
BOOST_AUTO_TEST_CASE(ThreeDMoments){
	const treecode::Node<3>& root = tree->getRoot();
	int total_charge, modulus_charge;
	V3d anal_coc  = anal.analytic_centre_of_charge(parts, total_charge, modulus_charge);
	V3d anal_dip  = anal.analytic_dipole_moment(parts, root.getCentreOfCharge());
	M3d anal_quad = anal.analytic_quadrupole_moment(parts, root.getCentreOfCharge());

	EIGEN_REQUIRE_CLOSE(anal_coc, root.getCentreOfCharge(), TOLERANCE);
	EIGEN_REQUIRE_CLOSE(anal_dip, root.getDipoleMoments(), TOLERANCE);
	EIGEN_REQUIRE_CLOSE(anal_quad, root.getQuadrupoleMoments(), TOLERANCE);
}

//Compare analytic forces, and check forces from tree
BOOST_AUTO_TEST_CASE(ThreeDForces){
	treecode::potentials::CoulombForceThreeD<3> pot(0.0, bounds);
	treecode::BarnesHutMAC<3> mac(0.0, bounds);
	treecode::Particle<3>* part = parts.front();

	V3d tree_force = force_on_particle(part, pot, treecode::potentials::quadrupole, mac);
	V3d anal_force = anal.analytic_force(*part, parts);
	double tree_pot = potential_at_particle(part, pot, treecode::potentials::quadrupole, mac);
	double anal_pot = anal.analytic_potential(*part, parts);

	//Should be identical with theta = 0
	EIGEN_REQUIRE_CLOSE(anal_force, tree_force, TOLERANCE);
	BOOST_REQUIRE_CLOSE(tree_pot, anal_pot, TOLERANCE);

	//Now try with theta = 0.5
	mac.setTheta(0.5);
	V3d mp_force = force_on_particle(part, pot, treecode::potentials::monopole, mac);
	double mp_pot = potential_at_particle(part, pot, treecode::potentials::monopole, mac);

	V3d dp_force = force_on_particle(part, pot, treecode::potentials::dipole, mac);
	double dp_pot = potential_at_particle(part, pot, treecode::potentials::dipole, mac);

	V3d qp_force = force_on_particle(part, pot, treecode::potentials::quadrupole, mac);
	double qp_pot = potential_at_particle(part, pot, treecode::potentials::quadrupole, mac);

	BOOST_TEST_MESSAGE("Is the following error acceptable? (theta = 0.5)");
	BOOST_TEST_MESSAGE("Analytic force:   (" << anal_force[0] << ", " << anal_force[1] << ", " << anal_force[2]<< ")");
	BOOST_TEST_MESSAGE("Monopole force:   (" << mp_force[0] << ", " << mp_force[1] << ", " << mp_force[2] << ")");
	BOOST_TEST_MESSAGE("Dipole force:     (" << dp_force[0] << ", " << dp_force[1] << ", " << dp_force[2]<< ")");
	BOOST_TEST_MESSAGE("Quadrupole force: (" << qp_force[0] << ", " << qp_force[1] << ", " << qp_force[2] << ")");
	BOOST_TEST_MESSAGE("Analytic potential:   " << anal_pot);
	BOOST_TEST_MESSAGE("Monopole potential:   " << mp_pot);
	BOOST_TEST_MESSAGE("Dipole potential:     " << dp_pot);
	BOOST_TEST_MESSAGE("Quadrupole potential: " << qp_pot);
}

BOOST_AUTO_TEST_SUITE_END();
