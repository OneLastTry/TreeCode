/**
 * @file
 * @brief A sample program using open boundaries in 3D.
 *
 *
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/foreach.hpp>
#include <algorithm>
#include <ctime>

#include <distributions/dists.h>

#include <potentials/CoulombForceEField.h>
#include <bounds/CountingPeriodicBounds.h>
#include <opt_parser/OptionParser.h>
#include <Particle.h>
#include <Node.h>
#include <Tree.h>
#include <TimeIntegrator.h>
#include <macs/BarnesHutMAC.h>
#include <bounds/OpenBoundary.h>
#include <pushers/LeapfrogPusher.h>
#include <bounds/FixedOpenBoundary.h>

#include <Eigen/Dense>
#include <vector>
#include <string>

using namespace std;
using namespace treecode;
using namespace treecode::potentials;

namespace po = boost::program_options;

typedef Eigen::Vector3d Vec;

Vec getForceOnParticle(Particle<3>* p,
		const Tree<3>& tree,
		const Potential<3>& potential,
		potentials::Precision precision,
		const AcceptanceCriterion<3>& mac) {
	typedef std::vector<Node<3>*> interaction_list;
	Vec force = Vec::Zero();

	interaction_list ilist;
	tree.getInteractionList(*p, ilist, mac);
	for (typename interaction_list::iterator it = ilist.begin(); it < ilist.end(); it++) {
		Node<3>* n = *it;
		force += potential.getForce(*p, *n, precision);
	}
	return force;
}

vector<Vec> getForces(const vector<Particle<3>*>& parts,
		const vector<int>& particle_indices,
		const Tree<3>& tree,
		const Potential<3>& potential,
		Precision prec,
		const AcceptanceCriterion<3>& mac){
	vector<Vec> forces;

	for(unsigned int i = 0; i < particle_indices.size(); i++){
		int part_index = particle_indices[i];
		Particle<3>* p = parts[part_index];
		Vec force = getForceOnParticle(p, tree, potential, prec, mac);
		forces.push_back(force);
	}
	return forces;
}

double getForceError(const vector<Vec>& direct_forces, const vector<Vec>& approx_forces){
	Vec error = Vec::Zero();
	Vec force_squared = Vec::Zero();

	for(unsigned int i = 0; i < direct_forces.size(); i++){
		for(unsigned int c = 0; c < error.rows(); c++){
			error[c] += (approx_forces[i][c] - direct_forces[i][c])*(approx_forces[i][c] - direct_forces[i][c]);
			force_squared[c] += approx_forces[i][c]*approx_forces[i][c];
		}
	}

	double total_error = 0;
	for(int c=0;c<error.rows();c++){
		error[c] = sqrt((double)(error[c] / force_squared[c]));
		total_error += error[c] / error.rows();
	}
	return total_error;
}

void printTimings(){
	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;
	using boost::mt19937;

	mt19937 rng;
	double length = 10;
	Vec origin(-length / 2, -length / 2, -length / 2);
	Vec max(length / 2, length / 2, length / 2);
	UniformDistribution<mt19937> position_dist(rng, origin, max);
	ConstDistribution velocity_dist(Vec::Zero());
	ConstantChargeDistribution electron_charges(-1);
	ConstantChargeDistribution ion_charges(1);


	clock_t start_time;
	for(unsigned int n = 100; n < 100000; n+=100){
		vector<Particle<3>*> parts;
		vector<Particle<3>*> ions = Particle<3>::generateParticles<mt19937>(n, 1837, rng, position_dist,
				velocity_dist, ion_charges);
		vector<Particle<3>*> electrons = Particle<3>::generateParticles<mt19937>(n, 1, rng, position_dist,
				velocity_dist, electron_charges);

		parts.insert(parts.end(), ions.begin(), ions.end());
		parts.insert(parts.end(), electrons.begin(), electrons.end());

		OpenBoundary<3> bounds;
		bounds.init(parts);
		CoulombForceThreeD<3> potential(0.0, bounds);
		LeapfrogPusher<3> push(0.01, bounds, potential);
		BarnesHutMAC<3> mac(0.0, bounds);
		Tree<3> tree(bounds, parts);
		tree.rebuild();


		vector<int> curr_pindices;
		for(unsigned int i=0;i<n;i++)
			curr_pindices.push_back(i);
		cerr << curr_pindices.size() << "\t";

		start_time = clock();
		getForces(parts, curr_pindices, tree, potential, Precision::monopole, mac);
		cerr << (clock() - start_time) << "\t";

		start_time = clock();
		getForces(parts, curr_pindices, tree, potential, Precision::dipole, mac);
		cerr << (clock() - start_time) << "\t";

		start_time = clock();
		getForces(parts, curr_pindices, tree, potential, Precision::quadrupole, mac);
		cerr << (clock() - start_time) << endl;

		Particle<3>::deleteParticles(parts);
	}

}

void printErrors(unsigned int n_test, vector<Particle<3>*>& parts,
		const Tree<3>& tree,
		const Potential<3>& potential,
		BarnesHutMAC<3>& mac){
	//Generate the list of particles to look at
	//This is a list of random indices so that we can take
	//a random sample of the particles.
	boost::mt19937 rng;
	boost::random::uniform_01<double> dist;
	vector<int> particle_indices;

	//Add particles until we have n_test
	while(particle_indices.size() < n_test){
		int x = (int)(dist(rng) * parts.size());
		//Make sure it doesn't already contain this index.
		if(std::find(particle_indices.begin(), particle_indices.end(), x) == particle_indices.end())
			particle_indices.push_back(x);
	}


	clock_t start_time, end_time;

	//Generate the direct force sum
	start_time = clock();
	vector<Vec> direct_forces = getForces(parts, particle_indices, tree, potential, Precision::monopole, mac);
	end_time = clock();
	//Output data for start.
	cout << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t";
	cout << (end_time-start_time) << "\t" << (end_time-start_time) << "\t" << (end_time-start_time) << "\t" << endl;

	for(double theta = 0.01; theta < 2.0; theta += 0.01){
		mac.setTheta(theta);

		start_time = clock();
		vector<Vec> approx_forces = getForces(parts, particle_indices, tree, potential, Precision::monopole, mac);
		clock_t monopole_time = clock() - start_time;
		double monopole_error = getForceError(direct_forces, approx_forces);

		start_time = clock();
		approx_forces = getForces(parts, particle_indices, tree, potential, Precision::dipole, mac);
		clock_t dipole_time = clock() - start_time;
		double dipole_error = getForceError(direct_forces, approx_forces);

		start_time = clock();
		approx_forces = getForces(parts, particle_indices, tree, potential, Precision::quadrupole, mac);
		double quadrupole_error = getForceError(direct_forces, approx_forces);
		clock_t quadrupole_time = clock() - start_time;

		cout << theta << "\t" << monopole_error << "\t" << dipole_error << "\t" << quadrupole_error << "\t";
		cout << monopole_time << "\t" << dipole_time << "\t" << quadrupole_time << endl;
	}
}




//Estimate error in leapfrog pusher by doing 20 orbits, and finding difference at start and end
double getLeapfrogError(double timestep){
	using namespace treecode;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;
	using std::vector;
	typedef Eigen::Matrix<double, 2, 1> Vec;

	Eigen::Vector3d position(1,2,3);
	Eigen::Vector3d velocity(1,2,3);

	Vec initial_position(100,0);
	Particle<2> p1(1e6, 1e6, Vec::Zero(), Vec(0, 0));
	Particle<2> p2(-1, 1, initial_position, Vec(0, -100));
	vector<Particle<2>*> all_parts;
	vector<Particle<2>*> push_parts;
	all_parts.push_back(&p1);
	all_parts.push_back(&p2);
	push_parts.push_back(&p2);

	FixedOpenBoundary<2> bounds(Vec(-110, -110), 120);
	CoulombForceThreeD<2> potential(0.0, bounds);
	Tree<2> tree(bounds, all_parts);
	LeapfrogPusher<2> push(timestep, bounds, potential);
	BarnesHutMAC<2> mac(0.0, bounds);
	push.init(push_parts, tree, monopole, mac);
	TimeIntegrator<2> integrator(timestep, M_PI*20, push_parts, tree, bounds, push, mac);
	integrator.start(monopole, 1);

	double error = (p2.getPosition() - initial_position).norm() / p2.getPosition().norm();
	return error;
}

void printLeapfrogError(){
	std::ofstream lerrors("leapfrog_errors.csv");
	for(double dt = 0.0001; dt < 1; dt += 0.0001){
		lerrors << dt << "\t" << getLeapfrogError(dt) << std::endl;
	}
}

int main(int argc, char **argv) {
	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;

	using boost::mt19937;
	mt19937 rng;

	double length = 10;
	unsigned int num_particles = 1000;
	double force_softening = 0.0;
	double timestep = 0.01;

	Vec origin(-length / 2, -length / 2, -length / 2);
	Vec max(length / 2, length / 2, length / 2);
	UniformDistribution<mt19937> position_dist(rng, origin, max);
	ConstDistribution velocity_dist(Vec::Zero());
	ConstantChargeDistribution electron_charges(-1);
	ConstantChargeDistribution ion_charges(1);

	vector<Particle<3>*> parts;
	vector<Particle<3>*> ions = Particle<3>::generateParticles<mt19937>(num_particles, 1837, rng, position_dist,
			velocity_dist, ion_charges);
	vector<Particle<3>* > electrons = Particle<3>::generateParticles<mt19937>(num_particles, 1, rng, position_dist,
			velocity_dist, electron_charges);

	parts.insert(parts.end(), ions.begin(), ions.end());
	parts.insert(parts.end(), electrons.begin(), electrons.end());

	OpenBoundary<3> bounds;
	bounds.init(parts);
	CoulombForceThreeD<3> potential(force_softening, bounds);
	LeapfrogPusher<3> push(timestep, bounds, potential);
	BarnesHutMAC<3> mac(0.0, bounds);
	Tree<3> tree(bounds, parts);
	tree.rebuild();

	printLeapfrogError();
	printErrors(1000, parts, tree, potential, mac);
	printTimings();
	return 0;
}
