/*
 * IOTest.cpp
 *
 *  Created on: 30 Mar 2012
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
#include <io/CoordTracker.h>
#include <io/ParticleReader.h>

#ifndef TOLERANCE
#define TOLERANCE 1e-10
#endif

#include "custom_asserts.h"

/**
 * Simple fixture that just generates some particles.
 */
template <int D>
struct IOFixture {
	IOFixture():pos_file("ion_pos.csv"), vel_file("ion_vel.csv"){
		using boost::random::mt19937;
		using namespace treecode;
		using namespace treecode::distribution;
		typedef Eigen::Vector3d V;

		double length = 1.0;
		unsigned int num_particles = 1000;
		mt19937 rng;

		//Create distributions
		UniformDistribution<mt19937> 	position_dist(rng, V::Zero(), V::Zero().array() + length);
		ConstDistribution				velocity_dist(V::Zero());
		ConstantChargeDistribution		electron_charges(-1);
		ConstantChargeDistribution		ion_charges(1);

		//Generate some particles
		ions = treecode::Particle<D>::template generateParticles<mt19937>(num_particles, 1837, rng,
				position_dist, velocity_dist, ion_charges);
	}

	void write_particles(){
		using namespace treecode;
		output::CoordTracker<3> ion_pos_tracker(pos_file, ions, output::CoordTracker<3>::POSITION);
		output::CoordTracker<3> ion_vel_tracker(vel_file, ions, output::CoordTracker<3>::VELOCITY);
		ion_pos_tracker.output();
		ion_vel_tracker.output();
	}

	~IOFixture(){
		treecode::Particle<D>::deleteParticles(ions);
		std::remove(pos_file);
		std::remove(vel_file);
	}

	const char *pos_file;
	const char *vel_file;
	std::vector<treecode::Particle<D>* > ions;
};

BOOST_FIXTURE_TEST_SUITE(IOTest, IOFixture<3>)

//Write particles to file, then make sure they are read back correctly
BOOST_AUTO_TEST_CASE(WriteRead){
	using namespace treecode;
	write_particles();

	io::ParticleReader<3> ion_reader(pos_file, vel_file, 1837, 1);
	std::vector<Particle<3>* > read_parts = ion_reader.readParticles(0);

	//Check we read back the same number of particles
	BOOST_REQUIRE_EQUAL(read_parts.size(), ions.size());
	typedef Particle<3> part_t;
	//Check that the positions and velocities are the same
	for(int i=0;i<read_parts.size();i++){
		part_t *orig_p = ions[i];
		part_t *read_p = read_parts[i];
		EIGEN_REQUIRE_CLOSE(orig_p->getPosition(), read_p->getPosition(), TOLERANCE);
		EIGEN_REQUIRE_CLOSE(orig_p->getVelocity(), read_p->getVelocity(), TOLERANCE);
	}
}

//Test reading failures
BOOST_AUTO_TEST_CASE(IOFail){
	using namespace treecode;
	//Write some particles
	write_particles();

	//Require an exception when attempting to read a nonexistent file
	BOOST_REQUIRE_THROW(
			io::ParticleReader<3> ion_reader("nonexistent", vel_file, 1,1),
			io::ReadError
			);
	BOOST_REQUIRE_THROW(
			io::ParticleReader<3> ion_reader(pos_file, "nonexistent", 1,1),
			io::ReadError
			);

	//Write the velocity file again, with one fewer particle :
	std::vector<Particle<3>* > one_fewer;
	one_fewer.insert(one_fewer.end(), ions.begin(), ions.end()-1);
	{
		output::CoordTracker<3> ion_vel_tracker(vel_file, one_fewer, output::CoordTracker<3>::VELOCITY);
		ion_vel_tracker.output();
	}
	//Read it back
	io::ParticleReader<3> ion_reader(pos_file, vel_file, 1837, 1);
	BOOST_REQUIRE_THROW(
			ion_reader.readParticles(),
			io::ReadError
			)
	//Do it again, with one fewer position
	{
		write_particles();
		output::CoordTracker<3> ion_vel_tracker(vel_file, one_fewer, output::CoordTracker<3>::POSITION);
		ion_vel_tracker.output();
	}
	io::ParticleReader<3> ion_reader2(pos_file, vel_file, 1837, 1);
	BOOST_REQUIRE_THROW(
			ion_reader.readParticles(),
			io::ReadError
			);
}

BOOST_AUTO_TEST_SUITE_END();
