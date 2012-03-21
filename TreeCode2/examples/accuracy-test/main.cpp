/**
 * @file
 * @brief A sample program using open boundaries in 3D.
 *
 *
 */

#include <cstdlib>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/foreach.hpp>
#include <algorithm>
#include <ctime>

#include <3d_typedefs.h>

#include <potentials/CoulombForceEField.h>
#include <bounds/CountingPeriodicBounds.h>
#include <opt_parser/OptionParser.h>
#include <Configuration.h>
#include <Particle.h>
#include <Node.h>
#include <Tree.h>
#include <TimeIntegrator.h>
#include <macs/BarnesHutMAC.h>

#include <Eigen/Dense>
#include <vector>
#include <string>

using namespace std;
using namespace treecode;
using namespace treecode::potentials;

namespace po = boost::program_options;

Vec getForceOnParticle(Particle3d* p, const Tree<Vec, Mat>& tree,
const Potential<Vec, Mat>& potential,
potentials::Precision precision, const AcceptanceCriterion<Vec, Mat>& mac) {
	typedef Node<Vec, Mat> Node;
	typedef std::vector<Node*> interaction_list;
	Vec force = Vec::Zero();

	interaction_list ilist;
	tree.getInteractionList(*p, ilist, mac);
	for (typename interaction_list::iterator it = ilist.begin(); it < ilist.end(); it++) {
		Node* n = *it;
		force += potential.getForce(*p, *n, precision);
	}
	return force;
}

vector<Vec> getForces(const vector<Particle3d*>& parts, const vector<int>& particle_indices, const Tree3d& tree, const Potential<Vec, Mat>& potential, Precision prec, const AcceptanceCriterion<Vec, Mat>& mac){
	vector<Vec> forces;

	for(unsigned int i = 0; i < particle_indices.size(); i++){
		int part_index = particle_indices[i];
		Particle3d* p = parts[part_index];
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
		error[c] = sqrt(error[c] / force_squared[c]);
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
	UniformDistribution3d position_dist(origin, max);
	ConstDistribution3d velocity_dist(Vec::Zero());
	ConstantChargeDistribution3d electron_charges(-1);
	ConstantChargeDistribution3d ion_charges(1);


	clock_t start_time;
	for(unsigned int n = 100; n < 100000; n+=100){
		int id = 0;
		vector<Particle3d*> parts;
		vector<Particle3d*> ions = Particle3d::generateParticles<mt19937>(n, 1837, rng, position_dist,
				velocity_dist, ion_charges, id);
		vector<Particle3d*> electrons = Particle3d::generateParticles<mt19937>(n, 1, rng, position_dist,
				velocity_dist, electron_charges, id);

		parts.insert(parts.end(), ions.begin(), ions.end());
		parts.insert(parts.end(), electrons.begin(), electrons.end());

		OpenBoundary3d bounds;
		bounds.init(parts);
		CoulombForce3d potential(0.0, bounds);
		LeapfrogPusher3d push(0.01, bounds, potential);
		BarnesHutMAC<Vec, Mat> mac(0.0, bounds);
		Tree3d tree(bounds, parts);
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

		Particle3d::deleteParticles(parts);
	}

}

void printErrors(unsigned int n_test, vector<Particle3d*>& parts, const Tree3d& tree, const Potential<Vec,Mat>& potential, BarnesHutMAC<Vec, Mat>& mac){
	//Generate the list of particles to look at
	boost::mt19937 rng;
	boost::uniform_01<double> dist;
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
	UniformDistribution3d position_dist(origin, max);
	ConstDistribution3d velocity_dist(Vec::Zero());
	ConstantChargeDistribution3d electron_charges(-1);
	ConstantChargeDistribution3d ion_charges(1);

	int id = 0;
	vector<Particle3d*> parts;
	vector<Particle3d*> ions = Particle3d::generateParticles<mt19937>(num_particles, 1837, rng, position_dist,
			velocity_dist, ion_charges, id);
	vector<Particle3d*> electrons = Particle3d::generateParticles<mt19937>(num_particles, 1, rng, position_dist,
			velocity_dist, electron_charges, id);

	parts.insert(parts.end(), ions.begin(), ions.end());
	parts.insert(parts.end(), electrons.begin(), electrons.end());

	OpenBoundary3d bounds;
	bounds.init(parts);
	CoulombForce3d potential(force_softening, bounds);
	LeapfrogPusher3d push(timestep, bounds, potential);
	BarnesHutMAC<Vec,Mat> mac(0.0, bounds);
	Tree3d tree(bounds, parts);
	tree.rebuild();

	printErrors(1000, parts, tree, potential, mac);
	printTimings();
	return 0;
}
