/*
 * Simulator.cpp
 *
 *  Created on: 30 Mar 2012
 *      Author: stefans
 */


#include <OptionParser.h>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <io/ParticleReader.h>
#include <boost/foreach.hpp>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

#include <io/ParticleReader.h>

//Boundary conditions
#include <bounds/OpenBoundary.h>
#include <bounds/PeriodicBoundary.h>
//Pusher
#include <pushers/LeapfrogPusher.h>
//Potentials
#include <potentials/CoulombForce.h>
#include <potentials/EwaldForce.h>
#include <potentials/InterpolatedEwaldSum.h>
#include <potentials/CoulombForceEField.h>
#include <potentials/DampingCoulombForce.h>
//MAC
#include <macs/BarnesHutMAC.h>
//Integrator
#include <TimeIntegrator.h>
//Output
#include <io/CoordTracker.h>
//General
#include <Particle.h>
#include <Tree.h>


#include <boost/preprocessor/iteration/local.hpp>

#define READ_PARTICLES_ND(z, n, unused) \
	int WOOPWOOPWOOP_##n(){\
		return n;\
	}

#define CALL_FUNCTIONS(z, n, unused) \
	DOSOMETHING_##n##d();

#define BOOST_PP_LOCAL_MACRO(n)   READ_PARTICLES_ND(~, n, ~)
#define BOOST_PP_LOCAL_LIMITS     (0, 1)
#include BOOST_PP_LOCAL_ITERATE()


namespace po = boost::program_options;

typedef std::vector<std::string> stringlist;
typedef std::vector<double> doublelist;
typedef std::vector<int>	intlist;

/**
 * Output generic parameters
 * @param theta
 * @param dt
 * @param max_time
 * @param force_softening
 */
void output_common_params(double theta, double dt, double max_time, double force_softening, const std::string& energy){
	std::cout << std::left << std::setw(20) << "Theta:" << theta << std::endl;
	std::cout << std::left << std::setw(20) << "Timestep:" << dt << std::endl;
	std::cout << std::left << std::setw(20) << "Max time:" << max_time << std::endl;
	std::cout << std::left << std::setw(20) << "Force softening:" << force_softening << std::endl;
	std::cout << std::left << std::setw(20) << "Energy file:" << energy << std::endl;
}

/**
 * @brief Read input files and bind output files.
 * @param opts	Options
 * @param part_it		Particles will be added to this.
 * @param tracker_it	Particle trackers will be added to this.
 */
typedef std::vector<treecode::Particle<3>* > part_list_3d;
typedef std::vector<treecode::output::CoordTracker<3>* > tracker_list_3d;
void read_particles_3d(const OptionParser& opts, part_list_3d& parts, tracker_list_3d& trackers){
	using treecode::output::CoordTracker;

	//Get the file lists, charges and masses
	stringlist pos_files = opts.get<stringlist>("pos-files");
	stringlist vel_files = opts.get<stringlist>("vel-files");
	doublelist masses = opts.get<doublelist>("masses");
	intlist charges = opts.get<intlist>("charges");
	//List of output files
	stringlist pos_out = opts.get<stringlist>("pos-out");
	stringlist vel_out = opts.get<stringlist>("vel-out");

	//Check vectors are all the same size
	int len = pos_files.size();
	if(vel_files.size() != len
			|| masses.size() != len
			|| charges.size() != len
			|| pos_out.size() != len
			|| vel_out.size() != len){
		std::cerr << "Must supply the same number of position files, velocity files, charges and masses." << std::endl;
		exit(1);
	}

	//Loop over each input file, generating particles
	for(int i=0;i<len;i++){
		treecode::io::ParticleReader<3> preader(pos_files[i], vel_files[i], masses[i], charges[i]);
		//Read the particles in that file
		part_list_3d current_parts = preader.readParticles();
		//Add them to the list of all particles
		parts.insert(parts.end(), current_parts.begin(), current_parts.end());
		if(opts.count("verbose")){
			std::cout << "Read " << current_parts.size() << " particles with m="
			<< masses[i] << ", q="
			<< charges[i] << " from "
			<< pos_files[i] << " and "
			<< vel_files[i] << std::endl;
		}

		//Bind them to a CoordTracker
		if(pos_out[i] != "none")
			trackers.push_back(new CoordTracker<3>(pos_out[i], current_parts, CoordTracker<3>::POSITION));
		if(vel_out[i] != "none")
			trackers.push_back(new CoordTracker<3>(vel_out[i], current_parts, CoordTracker<3>::VELOCITY));
	}
}

void simulate_open_3d(const OptionParser& opts){
	using namespace treecode;
	part_list_3d parts;
	tracker_list_3d trackers;

	read_particles_3d(opts, parts, trackers);

	//Read out some params we'll need
	double theta = opts.get<double>("theta");
	double force_softening = opts.get<double>("force-softening");
	double timestep = opts.get<double>("timestep");
	double max_time = opts.get<double>("max-time");
	int output_every = opts.get<int>("output-every");
	std::string energy_file = opts.get<std::string>("energy-out");

	if(opts.count("verbose")){
		std::cout << "Running a simulation with open boundaries" << std::endl;
		output_common_params(theta, timestep, max_time, force_softening, energy_file);
	}

	//Set up simulation
	OpenBoundary<3>		bounds;
	bounds.init(parts);

	BarnesHutMAC<3>						mac(theta, bounds);
	potentials::CoulombForceThreeD<3> 	*potential;
	if(opts.count("electric-field"))
		potential = new potentials::CoulombForceEField<3>(force_softening, bounds, opts.get<Eigen::VectorXd>("electric-field"));
	else if(opts.count("damping"))
		potential = new potentials::DampingCoulombForce<3>(force_softening, bounds, opts.get<double>("damping"));
	else
		potential = new potentials::CoulombForceThreeD<3>(force_softening, bounds);

	pusher::LeapfrogPusher<3> 			pusher(timestep, bounds, *potential);
	Tree<3>								tree(bounds, parts);
	TimeIntegrator<3>					integrator(timestep, max_time, parts, tree, bounds, pusher, mac);
	integrator.setEnergyOutputFile(energy_file.c_str());
	integrator.addParticleTrackers(trackers.begin(), trackers.end());

	if(opts.count("verbose"))
		std::cout << "Initialising pusher" << std::endl;
	pusher.init(parts, tree, potentials::quadrupole, mac);

	if(opts.count("verbose"))
		std::cout << "Starting time integrator" << std::endl;
	integrator.start(potentials::quadrupole, output_every);

	delete potential;
}

void simulate_periodic_3d(const OptionParser& opts){
	using namespace treecode;
	part_list_3d parts;
	tracker_list_3d trackers;
	read_particles_3d(opts, parts, trackers);

	//Read out some params we'll need
	double theta = opts.get<double>("theta");
	double force_softening = opts.get<double>("force-softening");
	double timestep = opts.get<double>("timestep");
	double max_time = opts.get<double>("max-time");
	Eigen::Vector3d origin = opts.get<Eigen::VectorXd>("origin");
	double length = opts.get<double>("length");
	int output_every = opts.get<int>("output-every");
	std::string energy_file = opts.get<std::string>("energy-out");

	if(opts.count("verbose")){
		std::cout << "Running a simulation with periodic boundaries" << std::endl;
		output_common_params(theta, timestep, max_time, force_softening, energy_file);
		std::cout << std::left << std::setw(20) << "Origin:";
		for(int i=0;i<origin.rows();i++)
			std::cout << origin[i] << " ";
		std::cout << std::endl;
		std::cout << std::left << std::setw(20) << "Length:" << length << std::endl;
	}

	//Fourier space and real space iterations
	int rs_its = opts.get<int>("real-space-iterations");
	int fs_its = opts.get<int>("fourier-space-iterations");

	PeriodicBoundary<3>		bounds(origin, length);
	BarnesHutMAC<3>		mac(theta, bounds);
	potentials::CoulombForceThreeD<3> *open_pot;
	if(opts.count("electric-field"))
		open_pot = new potentials::CoulombForceEField<3>(force_softening, bounds, opts.get<Eigen::VectorXd>("electric-field"));
	else if(opts.count("damping"))
		open_pot = new potentials::DampingCoulombForce<3>(force_softening, bounds, opts.get<double>("damping"));
	else
		open_pot = new potentials::CoulombForceThreeD<3>(force_softening, bounds);


	potentials::EwaldForce<3>			periodic_pot(force_softening, bounds, 2.0 / length, rs_its, fs_its);
	potentials::InterpolatedEwaldSum<3>	potential(force_softening, bounds, 40, periodic_pot, *open_pot);
	if(opts.count("verbose"))
		std::cout << "Initialising potential grid" << std::endl;
	potential.init();

	pusher::LeapfrogPusher<3> 		push(timestep, bounds, potential);
	Tree<3>					tree(bounds, parts);
	TimeIntegrator<3>		integrator(timestep, max_time, parts, tree, bounds, push, mac);
	integrator.setEnergyOutputFile(energy_file.c_str());
	integrator.addParticleTrackers(trackers.begin(), trackers.end());

	if(opts.count("verbose"))
		std::cout << "Initialising pusher" << std::endl;
	push.init(parts, tree, potentials::quadrupole, mac);

	typedef treecode::output::CoordTracker<3> tracker_t;
	BOOST_FOREACH(tracker_t* tr, trackers){
		tr->output();
	}

	if(opts.count("verbose"))
		std::cout << "Starting time integrator" << std::endl;
	integrator.start(potentials::quadrupole, output_every);

	delete open_pot;
}

int main(int argc, char **argv){
	using std::cout;
	using std::endl;
	Eigen::VectorXd origin;

	stringlist pos_file, vel_file;
	doublelist masses;
	intlist charges;

	OptionParser opts("Options");
	po::options_description gen_opts("Generic options");
	gen_opts.add_options()
		("verbose", 	"Increase verbosity of output.")
		("dimensions,d", po::value<int>()->required(), 			"Number of dimensions (2 or 3).")
		("timestep,t",	po::value<double>()->required(),		"Timestep (delta t).")
		("max-time,T",	po::value<double>()->required(),		"Maximum time to run until.")
		("bounds,b",	po::value<std::string>()->required(), 	"'open' or 'periodic'.")
		("theta,o",		po::value<double>()->required(),		"Critical opening angle.")
		("force-softening,f", po::value<double>()->required(), 	"Force softening constant.")
		("electric-field,E",po::value<Eigen::VectorXd>(), 		"Electric field (if any)")
		("damping,D",		po::value<double>(),				"Damping (if any)")
		("pos-files,p", 	po::value<stringlist>()->multitoken()->required(), 	"File containing positions of particles.")
		("vel-files,v", 	po::value<stringlist>()->multitoken()->required(), 	"File containing velocities of particles.")
		("pos-out",			po::value<stringlist>()->multitoken()->required(),	"Position output files")
		("vel-out",			po::value<stringlist>()->multitoken()->required(), 	"Velocity output files")
		("energy-out", 		po::value<std::string>()->required(), 				"Energy output file")
		("masses,m",	po::value<doublelist>()->multitoken()->required(), 	"List of particle masses.")
		("charges,q",	po::value<intlist>()->multitoken()->required(), 	"List of charges.")
		("output-every,e", po::value<int>()->default_value(1), 	"Output every arg timesteps")
	;

	po::options_description periodic_opts("Periodic options");
	periodic_opts.add_options()
		("origin,O", 	po::value<Eigen::VectorXd>(&origin), 	"Origin of the system.")
		("length,L",	po::value<double>(),					"Length of the system.")
		("real-space-iterations", po::value<int>()->default_value(5), "Number of Ewald iterations in real space")
		("fourier-space-iterations", po::value<int>()->default_value(5), "Number of Ewald iterations in fourier space")
	;
	opts.add(gen_opts);
	opts.add(periodic_opts);
	opts.setUsesConfigFile(true);
	opts.parse(argc, argv);

	std::string bounds_type = opts.get<std::string>("bounds");
	if(bounds_type == "open")
		simulate_open_3d(opts);
	else if(bounds_type == "periodic")
		simulate_periodic_3d(opts);
	else{
		std::cerr << "Bounds must be 'periodic' or 'open'." << std::endl;
		return 1;
	}

	return 0;
}
