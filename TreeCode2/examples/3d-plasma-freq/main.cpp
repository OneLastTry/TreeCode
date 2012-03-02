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

#include <Eigen/Dense>
#include <vector>
#include <string>

using namespace std;
using namespace treecode;

namespace po = boost::program_options;



Configuration3d parse_cmd_line(int argc, char **argv, double& length, unsigned int& num_parts, double& proportion, double& wavelengths, double& temperature){

	Option length_opt("--length", "-l", "Length of system (in debye lengths)");
	Option timestep_opt("--timestep", "-dt", "Individual timestep.");
	Option force_softening_opt("--softening", "-fs", "Force softening parameter.");
	Option max_time_opt("--max-time", "-mt", "Maximum time to run to.");
	Option theta_opt("--theta", "-t", "Critical opening angle.");
	Option param_opt("--param", "-p", "Plasma parameter (larger implies more ideal)");
	Option number_opt("--number", "-n", "Number of each species.");
	Option proportion_opt("--proportion", "-prop", "Proportion of particles in sinusoidal perturbation (0.0-1.0).");
	Option wavelength_opt("--wavelengths", "-w", "Number of wavelengths in perturbation.");
	Option temperature_opt("--temperature", "-temp", "Temperature of plasma.");
	Option help_opt("--help", "-h", "This help text.");
	Option opts[] = {
			length_opt, timestep_opt, force_softening_opt, max_time_opt, theta_opt, param_opt, number_opt,
			proportion_opt, wavelength_opt, temperature_opt, help_opt
	};

	OptionParser op(argc, argv);
	if(op.optionPresent(help_opt)){
		BOOST_FOREACH(Option o, opts) o.display(std::cout);
	}

	length = op.getOption<double>(length_opt);
	num_parts = op.getOption<double>(number_opt);
	proportion = op.getOption<double>(proportion_opt);
	wavelengths = op.getOption<double>(wavelength_opt);
	temperature = op.getOption<double>(temperature_opt);

	return Configuration3d(3,
			op.getOption<double>(theta_opt),
			op.getOption<double>(timestep_opt),
			op.getOption<double>(max_time_opt),
			op.getOption<unsigned int>(param_opt),
			op.getOption<double>(force_softening_opt));

}

int main(int argc, char **argv) {
	using namespace treecode::distribution;
	using namespace treecode::potentials;
	using namespace treecode::pusher;

	using boost::mt19937;
	mt19937 rng;


	double length, wavelengths, proportion, temperature;
	unsigned int num_particles;
	std::string dbname;
	Configuration3d c = parse_cmd_line(argc, argv, length, num_particles, proportion, wavelengths, temperature);

	Vec origin(-length/2,-length/2,-length/2);
	Vec max(length/2, length/2, length/2);
	UniformDistribution3d 			position_dist(origin, max);
	SinusoidalDistribution3d		perturbed_dist(0, origin, max, wavelengths, M_PI_2);
	ConstDistribution3d				i_velocity_dist(Vec::Zero());
	MaxwellDistribution3d			e_velocity_dist(1, temperature);
	ConstantChargeDistribution3d	electron_charges(-1);
	ConstantChargeDistribution3d	ion_charges(1);

	int id = 0;
	vector<Particle3d*> parts;
	vector<Particle3d*> ions = Particle3d::generateParticles<mt19937>(num_particles, 100000, rng,
			position_dist, i_velocity_dist, ion_charges, id);
	vector<Particle3d*> electrons = Particle3d::generateParticles<mt19937>(num_particles * (1 - proportion), 1, rng,
				position_dist, e_velocity_dist, electron_charges, id);
	vector<Particle3d*> perturbed_electrons = Particle3d::generateParticles<mt19937>(num_particles*proportion, 1, rng,
					perturbed_dist, e_velocity_dist, electron_charges, id);

	parts.insert(parts.end(), ions.begin(), ions.end());
	parts.insert(parts.end(), electrons.begin(), electrons.end());
	parts.insert(parts.end(), perturbed_electrons.begin(), perturbed_electrons.end());

	PeriodicBoundary3d		bounds(c, origin, length);
	CoulombForce3d 			open_pot(c, bounds);
	EwaldForce3d			periodic_pot(c, bounds, 2.0 / length, 5, 5);
	InterpolatedEwaldSum3d	potential(c, bounds, 25, periodic_pot, open_pot);
	potential.init();
	LeapfrogPusher3d 		push(c, bounds, potential);
	Tree3d					tree(c, bounds, parts);
	TimeIntegrator3d		integrator(c, parts, tree, bounds, push);
	integrator.setEnergyOutputFile("energies.csv");
	integrator.setPositionOutputFile("positions.csv");
	integrator.setVelocityOutputFile("velocities.csv");

	cout << "Initialising pusher" << endl;

	push.init(parts, tree, quadrupole);

	integrator.start(quadrupole, 1);

	return 0;
}
