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

#include <3d_typedefs.h>

#include <opt_parser/OptionParser.h>
#include <Configuration.h>
#include <Particle.h>
#include <Node.h>
#include <Tree.h>
#include <TimeIntegrator.h>
#include <macs/BarnesHutMAC.h>
#include <output/CoordTracker.h>

#include <Eigen/Dense>
#include <vector>
#include <string>

using namespace std;
using namespace treecode;

namespace po = boost::program_options;



void parse_cmd_line(int argc, char **argv,
		double& timestep, double &theta, double &max_time, double &fs,
		double& length, unsigned int& num_parts, double& proportion, double& wavelengths, double& temperature, bool& is_periodic){

	OptionParser op(argc, argv);
	bool help = false, periodic = false, open = false;

	op << new ArgOption<double>("--length", "-l", "Length of system (in debye lengths)", length) <<
			new ArgOption<double>("--timestep", "-dt", "Individual timestep.", timestep) <<
			new ArgOption<double>("--softening", "-fs", "Force softening parameter.", fs) <<
			new ArgOption<double>("--max-time", "-mt", "Maximum time to run to.", max_time) <<
			new ArgOption<double>("--theta", "-t", "Critical opening angle.", theta) <<
			new ArgOption<unsigned int>("--number", "-n", "Number of each species.", num_parts) <<
			new ArgOption<double>("--proportion", "-prop", "Proportion of particles in sinusoidal perturbation (0.0-1.0).", proportion) <<
			new ArgOption<double>("--wavelengths", "-w", "Number of wavelengths in perturbation.", wavelengths) <<
			new ArgOption<double>("--temperature", "-temp", "Temperature of plasma.", temperature) <<
			new BoolOption("--open", "-o", "Use open boundary conditions", open) <<
			new BoolOption("--periodic", "-p", "Use periodic boundary conditions", periodic) <<
			new BoolOption("--help", "-h", "This help text.", help);
	unsigned int num_options = op.parse();

	if(!periodic && !open){
		op.display(std::cerr);
		exit(1);
	}else if(periodic && open){
		std::cerr << "Must specify only one of --periodic or --open." << std::endl;
		exit(1);
	}
	is_periodic = periodic;

	if(help || num_options != op.size() - 2){
		op.display();
		exit(0);
	}
}

int main(int argc, char **argv) {
	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;
	using namespace treecode::output;

	using boost::mt19937;
	mt19937 rng;

	double length = 1;
	double force_softening = 0.001;
	double theta = 0.0;

	Vec origin(-length/2,-length/2,-length/2);
	Vec max(length/2, length/2, length/2);

	int id = 0;
	vector<Particle3d*> parts;
	Vec part_pos(-length/2, -length/2, 0);
	Particle3d system_particle(1, 1, part_pos, Vec::Zero(), id++);
	Particle3d test_particle(1, 1, Vec::Zero(), Vec::Zero(), id++);
	parts.push_back(&system_particle);
	parts.push_back(&test_particle);

	PeriodicBoundary3d		bounds(origin, length);
	BarnesHutMAC<D>		mac(theta, bounds);
	CoulombForce3d 			open_pot(force_softening, bounds);
	EwaldForce3d			periodic_pot(force_softening, bounds, 2.0 / length, 5, 5);
	InterpolatedEwaldSum3d	potential(force_softening, bounds, 20, periodic_pot, open_pot);
	potential.init();
	Tree3d					tree(bounds, parts);

	typedef Node<D> Node;
	double jitter = 0.005;
	int divs = 25;
	double test_max = (length - 2*jitter) - (length/2-jitter);
	for(double x = -length/2 + jitter; x < test_max; x+= test_max/divs){
		for(double y = -length/2 + jitter; y < test_max; y+= test_max/divs){
			Vec pos(x,y,0);
			test_particle.setPosition(pos);
			tree.rebuild();
			vector<Node* > ilist;
			tree.getInteractionList(test_particle, ilist, mac);
			double p = 0;
			BOOST_FOREACH(Node *n, ilist)
				p += potential.getPotential(test_particle, *n, monopole);
			cout << pos[0] << "\t" << pos[1] << "\t" << p << endl;
		}
		cout << endl;
	}

	return 0;
}
