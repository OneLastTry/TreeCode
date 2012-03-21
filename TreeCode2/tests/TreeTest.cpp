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
    	typedef Particle<V,M> Particle;
    	mt19937 rng;

    	double length = 1.0;
    	unsigned int num_particles = 1000;

    	UniformDistribution<mt19937,V> 			position_dist(V::Zero(), V::Zero().array() + length);
    	ConstDistribution<mt19937,V>			velocity_dist(V::Zero());
    	ConstantChargeDistribution<mt19937>		electron_charges(-1);
    	ConstantChargeDistribution<mt19937>		ion_charges(1);

    	int id = 0;
    	vector<Particle*> ions = Particle::template generateParticles<mt19937>(num_particles, 1837, rng,
    			position_dist, velocity_dist, ion_charges, id);
    	vector<Particle*> electrons = Particle::template generateParticles<mt19937>(num_particles, 1, rng,
    	    			position_dist, velocity_dist, electron_charges, id);
    	parts.insert(parts.end(), ions.begin(), ions.end());
    	parts.insert(parts.end(), electrons.begin(), electrons.end());
    	bounds.init(parts);
    	tree = new Tree<V,M>(bounds, parts);
    	tree->rebuild();
    }

    ~F(){
    	delete tree;
    }

    OpenBoundary<V,M> bounds;
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
		if(p == &test_part)
			continue;
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
		if(p == &test_part)
			continue;
		V r = (test_part.getPosition() - p->getPosition());
		potential += p->getCharge() / r.norm();
	}
	return potential;
}

template<typename V, typename M>
V force_on_particle(
		Particle<V,M>* p,
		const Tree<V,M>& tree,
		const potentials::Potential<V, M>& potential,
		potentials::Precision precision,
		const AcceptanceCriterion<V, M>& mac
		) {
	typedef Node<V, M> Node;
	typedef std::vector<Node*> interaction_list;
	V force = V::Zero();

	interaction_list ilist;
	tree.getInteractionList(*p, ilist, mac);
	for (typename interaction_list::iterator it = ilist.begin(); it < ilist.end(); it++) {
		Node* n = *it;
		force += potential.getForce(*p, *n, precision);
	}
	return force;
}

template<typename V, typename M>
double potenital_at_particle(
		Particle<V,M>* p,
		const Tree<V,M>& tree,
		const potentials::Potential<V, M>& potential,
		potentials::Precision precision,
		const AcceptanceCriterion<V, M>& mac
		) {
	typedef Node<V, M> Node;
	typedef std::vector<Node*> interaction_list;
	double pot = 0;

	interaction_list ilist;
	tree.getInteractionList(*p, ilist, mac);
	for (typename interaction_list::iterator it = ilist.begin(); it < ilist.end(); it++) {
		Node* n = *it;
		pot += potential.getPotential(*p, *n, precision);
	}
	return pot;
}

typedef Eigen::Vector2d V2d;
typedef Eigen::Matrix2d M2d;
typedef F<V2d, M2d> Fixture2d;
BOOST_FIXTURE_TEST_SUITE(Moments, Fixture2d)

BOOST_AUTO_TEST_CASE(TwoDMoments){
	const Node<V2d, M2d>& root = tree->getRoot();
	int total_charge, modulus_charge;
	V2d anal_coc  = analytic_centre_of_charge<V2d,M2d>(parts, total_charge, modulus_charge);
	V2d anal_dip  = analytic_dipole_moment<V2d,M2d>(parts, root.getCentreOfCharge());
	M2d anal_quad = analytic_quadrupole_moment<V2d,M2d>(parts, root.getCentreOfCharge());

	for(unsigned int i=0;i<anal_coc.rows();i++)
		BOOST_CHECK_CLOSE(anal_coc[i], root.getCentreOfCharge()[i], TOLERANCE);
	for(unsigned int i=0;i<anal_dip.rows();i++)
		BOOST_CHECK_CLOSE(anal_dip[i], root.getDipoleMoments()[i], TOLERANCE);
	for(unsigned int i=0;i<anal_quad.rows();i++)
		for(unsigned int j=0;j<anal_quad.cols();j++)
			BOOST_CHECK_CLOSE(anal_quad.coeff(i,j), root.getQuadrupoleMoments().coeff(i,j), TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TwoDForces){
	using namespace potentials;
	CoulombForceThreeD<V2d,M2d> pot(0.0, bounds);
	BarnesHutMAC<V2d,M2d> mac(0.0, bounds);
	Particle<V2d,M2d>* part = parts.front();

	V2d tree_force = force_on_particle<V2d,M2d>(part, *tree, pot, quadrupole, mac);
	V2d anal_force = analytic_force<V2d,M2d>(*part, parts);
	double tree_pot = potenital_at_particle<V2d,M2d>(part, *tree, pot, quadrupole, mac);
	double anal_pot = analytic_potential<V2d,M2d>(*part, parts);

	//Should be identical with theta = 0
	for(unsigned int i=0;i< tree_force.rows();i++)
		BOOST_CHECK_CLOSE(anal_force[i], tree_force[i], TOLERANCE);
	BOOST_CHECK_CLOSE(tree_pot, anal_pot, TOLERANCE);

	//Now try with theta = 0.5

}

BOOST_AUTO_TEST_SUITE_END();
