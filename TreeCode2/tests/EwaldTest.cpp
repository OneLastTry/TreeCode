/*
 * EwaldTest.cpp
 *
 *  Created on: 20 Jan 2012
 *      Author: stefans
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <vector>


#include "../3d_typedefs.h"

using namespace std;
using namespace treecode;
using namespace treecode::potentials;
using namespace treecode::distribution;
using namespace Eigen;

typedef Eigen::Vector3d V;
typedef Eigen::Matrix3d M;

BOOST_AUTO_TEST_SUITE(Ewald)

double direct_periodic_potential(const Particle3d& p1, const Particle3d& p2, int iterations, double L) {
	double potential = 0;
	for (int i = -iterations; i < iterations; i++) {
		for (int j = -iterations; j < iterations; j++) {
			for (int k = -iterations; k < iterations; k++) {
				Vector3d dv(i*L, j*L, k*L);
				dv -= (p2.getPosition() - p1.getPosition());
				potential += p1.getCharge() * p2.getCharge() / dv.norm();
			}
		}
	}
	return potential;
}

V direct_periodic_force(const Particle3d& p1, const Particle3d& p2, int iterations, double L) {
	Vector3d force = Vector3d::Zero();
	for (int i = -iterations; i < iterations; i++) {
		for (int j = -iterations; j < iterations; j++) {
			for (int k = -iterations; k < iterations; k++) {
				V dv(i*L, j*L, k*L);
				dv += (p1.getPosition() - p2.getPosition());
				force += p1.getCharge() * p2.getCharge() * dv / (dv.norm() * dv.squaredNorm());
			}
		}
	}
	return force;
}

//BOOST_AUTO_TEST_CASE(EwaldTest3dPotential) {
//	Configuration c(3, 0.0, 0.01, 1, 1, 0);
//	PeriodicBoundary bounds(c, Vector3d::Zero(), 1);
//	EwaldForce potential(bounds, 3, 5, 5);
//
//	Particle p1(c, +1, 1, Vector3d(0.4, 0.5, 0), Vector3d::Zero());
//	Particle p2(c, -1, 1, Vector3d(0.6, 0.5, 0), Vector3d::Zero());
//
//	Particle test_part(c, +1, 1, Vector3d(0,0,0), Vector3d::Zero());
//
//	for(double x = 0;x < bounds.getSize(); x+= 0.1) {
//		for(double y = 0; y < bounds.getSize(); y+= 0.1) {
//			test_part.setPosition(Vector3d(x,y,0));
//			double pot = 0;
//			pot += potential.real_space_pot(test_part, p1);
//			pot += potential.fourier_space_pot(test_part, p1);
//			pot += potential.real_space_pot(test_part, p2);
//			pot += potential.fourier_space_pot(test_part, p2);
//
//			double one_part_pot = p1.getCharge()*test_part.getCharge() / (p1.getPosition() - test_part.getPosition()).norm();
//			one_part_pot += p2.getCharge() *test_part.getCharge() / (p2.getPosition() - test_part.getPosition()).norm();
//
//			double direct_pot = direct_periodic_potential(p1, test_part, 10, bounds.getSize());
//			direct_pot += direct_periodic_potential(p2, test_part, 10, bounds.getSize());
//
//			cout << x << "\t" << y << "\t" << pot << "\t" << one_part_pot << "\t" << direct_pot << endl;
//		}
//		cout << endl;
//	}
//}

//BOOST_AUTO_TEST_CASE(EwaldTest3dForce){
//	Configuration c(3, 0.0, 0.01, 1, 1, 0);
//	PeriodicBoundary bounds(c, Vector3d::Zero(), 10);
//	EwaldForce potential(bounds, 0.01, 20, 0);
//
//	Particle p1(c, +1, 1, Vector3d(0, 0, 0), Vector3d::Zero());
//	Particle p2(c, -1, 1, Vector3d(0.6, 0.5, 0), Vector3d::Zero());
//
//	Particle test_part(c, +1, 1, Vector3d(0,0,0), Vector3d::Zero());
//
//	for(double x = 0.5;x < bounds.getSize() - 0.5; x+= 0.25) {
//		for(double y = 0.5; y < bounds.getSize() - 0.5; y+= 0.25) {
//			test_part.setPosition(Vector3d(x,y,0));
//			Vector3d force = Vector3d::Zero();
//			force += potential.real_space_force(test_part, p1);
//			force += potential.fourier_space_force(test_part, p1);
//
//			Vector3d one_part_force = Vector3d::Zero();
//			double r = (test_part.getPosition() - p1.getPosition()).norm();
//			double r_3 = r*r*r;
//			one_part_force += test_part.getCharge() * p1.getCharge() *
//					(test_part.getPosition() - p1.getPosition()) / r_3;
//
//			Vector3d direct_force = Vector3d::Zero();
//			direct_force += direct_periodic_force(test_part, p1, 10, bounds.getSize());
//
//			cout << x << "\t" << y << "\t" << force[0] << "\t" << force[1] << "\t" <<
//					one_part_force[0] << "\t" << one_part_force[1] << "\t" <<
//					direct_force[0] << "\t" << direct_force[1] << "\t" << endl;
//		}
//		cout << endl;
//	}
//}

//BOOST_AUTO_TEST_CASE(EwaldTestDipole){
//	using boost::random::mt19937;
//	int num_particles = 10;
//	mt19937 rng(0);
//
//	Configuration c(3, 0.0, 0.01, 1, 1, 0);
//	PeriodicBoundary bounds(c, Vector3d::Zero(), 50);
//	EwaldForce potential(bounds, 0.2, 10, 10);
//
//	ConstDistribution<mt19937> particle_velocities(Vector3d::Zero());
//	UniformDistribution<mt19937> particle_positions(c, Vector3d::Zero(), Vector3d(1, 1, 1));
//	ConstantChargeDistribution<mt19937> proton_charges(+1);
//	ConstantChargeDistribution<mt19937> electron_charges(-1);
//
//	vector<Particle*> parts;
//	vector<Particle*> protons = Particle::generateParticles<mt19937>(c, num_particles/2, 1837, rng,
//			particle_positions, particle_velocities, proton_charges);
//	vector<Particle*> electrons = Particle::generateParticles<mt19937>(c, num_particles/2, 1, rng,
//				particle_positions, particle_velocities, electron_charges);
//	parts.insert(parts.end(), protons.begin(), protons.end());
//	parts.insert(parts.end(), electrons.begin(), electrons.end());
//
//	Tree tr(c, bounds, parts);
//	tr.rebuild();
//
//	Particle test_particle(c, +1, 1, Vector3d(15,0,0), Vector3d(0,0,0));
//
////	double direct_pot = 0;
//	Vector3d ewald_force = Vector3d::Zero();
//	Vector3d direct_force = Vector3d::Zero();
//	for(Particle* p : parts){
//		ewald_force += potential.real_space_force(test_particle, *p);
//		ewald_force += potential.fourier_space_force(test_particle, *p);
//
//		direct_force += direct_periodic_force(test_particle, *p, 10, bounds.getSize());
////		direct_pot += potential.real_space_pot(test_particle, *p);
////		direct_pot += potential.fourier_space_pot(test_particle, *p);
//	}
//
////	double dipole_pot = potential.getPotential(test_particle, tr.getRoot(), Potential::Precision::dipole);
//	Vector3d dipole_force = potential.getForce(test_particle, tr.getRoot(), Potential::Precision::quadrupole);
//	cout << ewald_force[0] << "\t" << ewald_force[1] << "\t" << ewald_force[2] << "\t" << endl <<
//			direct_force[0] << "\t" << direct_force[1] << "\t" << direct_force[2] << "\t" << endl <<
//			dipole_force[0] << "\t" << dipole_force[1] << "\t" << dipole_force[2] << "\t" << endl;
////	cout << direct_pot << "\t" << dipole_pot << endl;
//}

//BOOST_AUTO_TEST_CASE(InterpolatedEwaldPotential){
//	Configuration c(3, 0.0, 0.01, 1, 1, 0);
//	PeriodicBoundary bounds(c, Vector3d::Zero(), 1);
//	EwaldForce pot(c, bounds, 2.0 / bounds.getSize(), 5, 5);
//
//	vector<Particle*> parts;
//	Particle* p1 = new Particle(c, 1, 1, Vector3d(0.025,0.025,0), Vector3d::Zero());
//	parts.push_back(p1);
//
//	Tree tr(c, bounds, parts);
//	tr.rebuild();
//
//	Particle* test_particle = new Particle(c, 1, 1, Vector3d::Zero(), Vector3d::Zero());
//
//	for(double x = 0; x < bounds.getSize()+0.05; x+=0.05){
//		for(double y = 0; y < bounds.getSize()+0.05; y+=0.05){
//			test_particle->setPosition(Vector3d(x,y,0));
//			cout << x << "\t" << y << "\t" << pot.getPotential(*test_particle, tr.getRoot(), Potential::quadrupole) << endl;
//		}
//		cout << endl;
//	}
//}

//BOOST_AUTO_TEST_CASE(InterpolatedEwaldPotential){
//	Configuration3d c(3, 0.0, 0.01, 1, 1, 0.1);
//	PeriodicBoundary3d bounds(c, V::Zero(), 5);
//
//	EwaldForce3d ewald_pot(c, bounds, 2.0 / bounds.getSize(), 5, 5);
//	CoulombForce3d coulomb_force(c, bounds);
//
//	InterpolatedEwaldSum3d pot(c, bounds, 20, ewald_pot, coulomb_force);
//	pot.init();
//
//	vector<Particle3d*> parts;
//	Particle3d* p1 = new Particle3d(-1, 1, Vec(0.025,0.025,0), Vec::Zero(), 1);
//	parts.push_back(p1);
//
//	Tree3d tr(c, bounds, parts);
//	tr.rebuild();
//
//	Particle3d* test_particle = new Particle3d(1, 1, Vec::Zero(), Vec::Zero(), 2);
//
//	pot.outputField();
//
//	ofstream pot_out("potential.csv");
//	for(double x = 0; x < bounds.getSize()+0.25; x+=0.25){
//		for(double y = 0; y < bounds.getSize()+0.25; y+=.025){
//			test_particle->setPosition(Vector3d(x,y,0));
//			pot_out << x << "\t" << y << "\t" << pot.getPotential(*test_particle, tr.getRoot(), potentials::quadrupole) << endl;
//		}
//		pot_out << endl;
//	}
//	pot_out.close();
//
//	ofstream force_out("force.csv");
//	for(double x = 0; x < bounds.getSize()+0.25; x+=0.25){
//		for(double y = 0; y < bounds.getSize()+0.25; y+=0.25){
//			test_particle->setPosition(Vector3d(x,y,0));
//			Vec force = pot.getForce(*test_particle, tr.getRoot(), potentials::quadrupole);
//			force_out << x << "\t" << y << "\t" << force[0] << "\t" << force[1] << endl;
//		}
//		force_out << endl;
//	}
//	force_out.close();
//}

//BOOST_AUTO_TEST_CASE(InterpolatedEwaldPotential){
//	using boost::random::mt19937;
//	Configuration3d c(3, 0.0, 0.01, 1, 1, 0.1);
//	PeriodicBoundary3d bounds(c, Vector3d::Zero(), 10);
//	OpenBoundary3d open_bounds(c);
//	mt19937 rng(0);
//
//	ConstDistribution3d particle_velocities(Vec::Zero());
//	UniformDistribution3d particle_positions(Vec::Zero(), Vec(1,1,1));
//	ConstantChargeDistribution3d proton_charges(+1);
//	ConstantChargeDistribution3d electron_charges(-1);
//
//	vector<Particle3d*> parts;
//	int num_particles = 5;
//	int id = 0;
//
//	vector<Particle3d*> protons = Particle3d::generateParticles<mt19937>(num_particles, 1837, rng,
//			particle_positions, particle_velocities, proton_charges, id);
//	vector<Particle3d*> electrons = Particle3d::generateParticles<mt19937>(num_particles, 1, rng,
//				particle_positions, particle_velocities, electron_charges, id);
//
//	parts.insert(parts.end(), protons.begin(), protons.end());
//	parts.insert(parts.end(), electrons.begin(), electrons.end());
//	Particle3d* test_particle = new Particle3d(1, 1, V::Zero(), V::Zero(), ++id);
//	parts.push_back(test_particle);
//
//
//	Tree3d tr(c, bounds, parts);
//	tr.rebuild();
//
//	EwaldForce3d ewald_pot(c, bounds, 2.0 / bounds.getSize(), 5, 5);
//	CoulombForce3d coulomb_force(c, open_bounds);
//
//	InterpolatedEwaldSum3d pot(c, bounds, 20, ewald_pot, coulomb_force);
//	pot.init();
//
//	pot.outputField();
//
//	ofstream pot_out("multipole_potential.csv");
//	for(double x = 0; x < bounds.getSize()+0.25; x+=0.25){
//		for(double y = 0; y < bounds.getSize()+0.5; y+=0.25){
//			test_particle->setPosition(Vector3d(x,y,0));
//			tr.rebuild();
//
//			double p = 0;
//			std::vector<Node3d*> ilist;
//			tr.getInteractionList(*test_particle, ilist);
//			for(Node3d* n : ilist)
//				p += pot.getPotential(*test_particle, *n, potentials::quadrupole);
//
//			pot_out << x << "\t" << y << "\t" << p << endl;
//		}
//		pot_out << endl;
//	}
//	pot_out.close();
//
//	ofstream force_out("multipole_force.csv");
//	for(double x = 0; x < bounds.getSize()+0.25; x+=0.25){
//		for(double y = 0; y < bounds.getSize()+0.25; y+=0.25){
//			test_particle->setPosition(Vector3d(x,y,0));
//			tr.rebuild();
//
//			std::vector<Node3d*> ilist;
//			tr.getInteractionList(*test_particle, ilist);
//			Vector3d force = Vector3d::Zero();
//			for(Node3d* n : ilist)
//				force += pot.getForce(*test_particle, *n, potentials::quadrupole);
//
//			force_out << x << "\t" << y << "\t" << force[0] << "\t" << force[1] << endl;
//		}
//		force_out << endl;
//	}
//	force_out.close();
//}

BOOST_AUTO_TEST_SUITE_END()
